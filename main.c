#include <raylib.h>
#include <raymath.h>
#include <math.h>
#include <stdio.h>

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600
#define LINK_LENGTH 100.0f
#define CYAN (Color){ 0, 255, 255, 255 }

#define MAX_HISTORY 200

typedef struct {
    float theta1; // Angle of first link relative to base
    float theta2; // Angle of second link relative to first
    float theta3; // Angle of third link relative to second
} RobotArm;

typedef struct {
    float J[2][3];
    float JJT[2][2];
    float manipulability;
    float eigenvalues[2];
    Vector2 eigenvectors[2];
} JacobianInfo;

// Global history for plotting
float lambda1History[MAX_HISTORY] = {0};
float lambda2History[MAX_HISTORY] = {0};
int historyIndex = 0;

// Calculate Jacobian Determinant (Manipulability Index)
JacobianInfo CalculateManipulability(RobotArm arm) {
    float L1 = LINK_LENGTH;
    float L2 = LINK_LENGTH;
    float L3 = LINK_LENGTH;

    float t1 = arm.theta1;
    float t2 = arm.theta1 + arm.theta2;
    float t3 = arm.theta1 + arm.theta2 + arm.theta3;

    // Partial derivatives of x
    float dx_d1 = -L1*sinf(t1) - L2*sinf(t2) - L3*sinf(t3);
    float dx_d2 = -L2*sinf(t2) - L3*sinf(t3);
    float dx_d3 = -L3*sinf(t3);

    // Partial derivatives of y
    float dy_d1 = L1*cosf(t1) + L2*cosf(t2) + L3*cosf(t3);
    float dy_d2 = L2*cosf(t2) + L3*cosf(t3);
    float dy_d3 = L3*cosf(t3);

    // Jacobian Matrix J = [dx_d1 dx_d2 dx_d3]
    //                     [dy_d1 dy_d2 dy_d3]
    
    // Calculate JJ^T (2x2 matrix)
    // A = [a b]
    //     [c d]
    
    float j11 = dx_d1, j12 = dx_d2, j13 = dx_d3;
    float j21 = dy_d1, j22 = dy_d2, j23 = dy_d3;

    float a = j11*j11 + j12*j12 + j13*j13;
    float b = j11*j21 + j12*j22 + j13*j23;
    float c = b; // Symmetric
    float d = j21*j21 + j22*j22 + j23*j23;

    float detJJT = a*d - b*c;
    
    if (detJJT < 0) detJJT = 0; // Should not happen theoretically for real matrices but float precision
    
    JacobianInfo info;
    info.J[0][0] = j11; info.J[0][1] = j12; info.J[0][2] = j13;
    info.J[1][0] = j21; info.J[1][1] = j22; info.J[1][2] = j23;
    
    info.JJT[0][0] = a; info.JJT[0][1] = b;
    info.JJT[1][0] = c; info.JJT[1][1] = d;
    
    info.manipulability = sqrtf(detJJT);

    // Eigen Decomposition of 2x2 Matrix
    // Trace T = a + d
    // Determinant D = ad - bc
    // Eigenvalues = (T +/- sqrt(T^2 - 4D)) / 2
    
    float T = a + d;
    float D = a*d - b*c;
    float discriminant = T*T - 4*D;
    if (discriminant < 0) discriminant = 0;
    
    float lambda1 = (T + sqrtf(discriminant)) / 2.0f;
    float lambda2 = (T - sqrtf(discriminant)) / 2.0f;
    
    info.eigenvalues[0] = lambda1;
    info.eigenvalues[1] = lambda2;
    
    // Eigenvectors
    // For lambda1: (a - lambda1)x + by = 0
    // If b != 0, y = -(a - lambda1)/b * x. Let x = 1.
    // If b == 0, matrix is diagonal. Eigenvectors are (1,0) and (0,1).
    
    if (fabsf(b) > 1e-6) {
        info.eigenvectors[0] = (Vector2){ 1.0f, -(a - lambda1)/b };
        info.eigenvectors[1] = (Vector2){ 1.0f, -(a - lambda2)/b };
    } else {
        info.eigenvectors[0] = (Vector2){ 1.0f, 0.0f };
        info.eigenvectors[1] = (Vector2){ 0.0f, 1.0f };
    }
    
    // Normalize eigenvectors
    float len1 = sqrtf(info.eigenvectors[0].x*info.eigenvectors[0].x + info.eigenvectors[0].y*info.eigenvectors[0].y);
    float len2 = sqrtf(info.eigenvectors[1].x*info.eigenvectors[1].x + info.eigenvectors[1].y*info.eigenvectors[1].y);
    
    if (len1 > 0) { info.eigenvectors[0].x /= len1; info.eigenvectors[0].y /= len1; }
    if (len2 > 0) { info.eigenvectors[1].x /= len2; info.eigenvectors[1].y /= len2; }
    
    return info;
}

int main(void)
{
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "KSV: Kinematic Singularity Visualizer");
    SetTargetFPS(60);

    RobotArm arm = { 0.0f, 0.0f, 0.0f };
    Vector2 basePos = { SCREEN_WIDTH / 2.0f, SCREEN_HEIGHT / 2.0f };

    while (!WindowShouldClose()) {
        // Update
        float speed = 2.0f * GetFrameTime();
        
        if (IsKeyDown(KEY_Q)) arm.theta1 += speed;
        if (IsKeyDown(KEY_A)) arm.theta1 -= speed;
        
        if (IsKeyDown(KEY_W)) arm.theta2 += speed;
        if (IsKeyDown(KEY_S)) arm.theta2 -= speed;
        
        if (IsKeyDown(KEY_E)) arm.theta3 += speed;
        if (IsKeyDown(KEY_D)) arm.theta3 -= speed;

        JacobianInfo info = CalculateManipulability(arm);
        float manipulability = info.manipulability;

        // Update History
        lambda1History[historyIndex] = info.eigenvalues[0];
        lambda2History[historyIndex] = info.eigenvalues[1];
        historyIndex = (historyIndex + 1) % MAX_HISTORY;

        // Draw
        BeginDrawing();
        ClearBackground(RAYWHITE);

        // Calculate joint positions
        Vector2 p0 = basePos;
        Vector2 p1 = { p0.x + LINK_LENGTH * cosf(arm.theta1), p0.y + LINK_LENGTH * sinf(arm.theta1) };
        Vector2 p2 = { p1.x + LINK_LENGTH * cosf(arm.theta1 + arm.theta2), p1.y + LINK_LENGTH * sinf(arm.theta1 + arm.theta2) };
        Vector2 p3 = { p2.x + LINK_LENGTH * cosf(arm.theta1 + arm.theta2 + arm.theta3), p2.y + LINK_LENGTH * sinf(arm.theta1 + arm.theta2 + arm.theta3) };

        // Determine color based on manipulability
        Color armColor = GREEN;
        if (manipulability < 10000.0f) armColor = YELLOW; // Adjusted thresholds
        if (manipulability < 2000.0f) armColor = RED;

        // Draw Arm
        DrawLineEx(p0, p1, 5.0f, armColor);
        DrawLineEx(p1, p2, 5.0f, armColor);
        DrawLineEx(p2, p3, 5.0f, armColor);

        DrawCircleV(p0, 8.0f, BLACK);
        DrawCircleV(p1, 6.0f, BLACK);
        DrawCircleV(p2, 6.0f, BLACK);
        DrawCircleV(p3, 4.0f, RED);

        // Draw Manipulability Ellipsoid
        // Axis lengths are sqrt(eigenvalues)
        float axis1 = sqrtf(info.eigenvalues[0]);
        float axis2 = sqrtf(info.eigenvalues[1]);
        
        // Scale factor for visualization (arbitrary to make it visible)
        float scale = 0.5f; 
        
        Vector2 v1 = info.eigenvectors[0];
        Vector2 v2 = info.eigenvectors[1];
        
        // Draw Principal Axes
        DrawLineEx(p3, (Vector2){ p3.x + v1.x * axis1 * scale, p3.y + v1.y * axis1 * scale }, 3.0f, MAGENTA);
        DrawLineEx(p3, (Vector2){ p3.x + v2.x * axis2 * scale, p3.y + v2.y * axis2 * scale }, 3.0f, CYAN);
        
        // Draw Ellipse (Approximation with lines)
        int segments = 32;
        float angleStep = 2.0f * PI / segments;
        Vector2 prevPoint;
        
        // Parametric equation: P(t) = center + (v1 * axis1 * cos(t) + v2 * axis2 * sin(t)) * scale
        for (int i = 0; i <= segments; i++) {
            float t = i * angleStep;
            float c = cosf(t);
            float s = sinf(t);
            
            Vector2 point = {
                p3.x + (v1.x * axis1 * c + v2.x * axis2 * s) * scale,
                p3.y + (v1.y * axis1 * c + v2.y * axis2 * s) * scale
            };
            
            if (i > 0) {
                DrawLineEx(prevPoint, point, 2.0f, Fade(ORANGE, 0.6f));
            }
            prevPoint = point;
        }

        // Draw UI
        DrawText("Controls:", 10, 10, 20, DARKGRAY);
        DrawText("Joint 1: Q/A", 10, 35, 20, DARKGRAY);
        DrawText("Joint 2: W/S", 10, 60, 20, DARKGRAY);
        DrawText("Joint 3: E/D", 10, 85, 20, DARKGRAY);

        char buffer[64];
        sprintf(buffer, "Manipulability Index (w): %.2f", manipulability);
        DrawText(buffer, 10, SCREEN_HEIGHT - 40, 20, BLACK);

        // Draw Health Bar
        float maxW = 50000.0f; // Approximate max manipulability for scaling
        float barWidth = (manipulability / maxW) * 200.0f;
        if (barWidth > 200.0f) barWidth = 200.0f;
        
        DrawRectangle(10, SCREEN_HEIGHT - 70, 200, 20, LIGHTGRAY);
        DrawRectangle(10, SCREEN_HEIGHT - 70, (int)barWidth, 20, armColor);
        DrawRectangleLines(10, SCREEN_HEIGHT - 70, 200, 20, DARKGRAY);

        // Draw Matrices
        int startX = SCREEN_WIDTH - 350;
        int startY = 10;
        
        DrawText("Jacobian Matrix J:", startX, startY, 20, DARKGRAY);
        char matBuf[256];
        sprintf(matBuf, "[ %6.1f %6.1f %6.1f ]\n[ %6.1f %6.1f %6.1f ]", 
            info.J[0][0], info.J[0][1], info.J[0][2],
            info.J[1][0], info.J[1][1], info.J[1][2]);
        DrawText(matBuf, startX, startY + 25, 20, BLACK);

        DrawText("JJ^T Matrix:", startX, startY + 80, 20, DARKGRAY);
        sprintf(matBuf, "[ %8.1f %8.1f ]\n[ %8.1f %8.1f ]", 
            info.JJT[0][0], info.JJT[0][1],
            info.JJT[1][0], info.JJT[1][1]);
        DrawText(matBuf, startX, startY + 105, 20, BLACK);

        DrawText("Eigenvalues (lambda):", startX, startY + 160, 20, DARKGRAY);
        sprintf(matBuf, "L1: %.1f  L2: %.1f", info.eigenvalues[0], info.eigenvalues[1]);
        DrawText(matBuf, startX, startY + 185, 20, BLACK);

        DrawText("Singular Values (sigma):", startX, startY + 215, 20, DARKGRAY);
        sprintf(matBuf, "S1: %.1f  S2: %.1f", sqrtf(info.eigenvalues[0]), sqrtf(info.eigenvalues[1]));
        DrawText(matBuf, startX, startY + 240, 20, BLACK);

        // Draw Graph
        int graphX = startX;
        int graphY = startY + 280;
        int graphW = 300;
        int graphH = 100;
        
        DrawRectangle(graphX, graphY, graphW, graphH, Fade(LIGHTGRAY, 0.3f));
        DrawRectangleLines(graphX, graphY, graphW, graphH, DARKGRAY);
        
        float maxVal = 50000.0f; // Scale for graph
        
        for (int i = 0; i < MAX_HISTORY - 1; i++) {
            int idx = (historyIndex + i) % MAX_HISTORY;
            int nextIdx = (historyIndex + i + 1) % MAX_HISTORY;
            
            float x1 = graphX + (float)i / MAX_HISTORY * graphW;
            float x2 = graphX + (float)(i + 1) / MAX_HISTORY * graphW;
            
            float y1_l1 = graphY + graphH - (lambda1History[idx] / maxVal) * graphH;
            float y2_l1 = graphY + graphH - (lambda1History[nextIdx] / maxVal) * graphH;
            
            float y1_l2 = graphY + graphH - (lambda2History[idx] / maxVal) * graphH;
            float y2_l2 = graphY + graphH - (lambda2History[nextIdx] / maxVal) * graphH;
            
            // Clamp values
            if (y1_l1 < graphY) y1_l1 = graphY; 
            if (y2_l1 < graphY) y2_l1 = graphY;
            
            if (y1_l2 < graphY) y1_l2 = graphY; 
            if (y2_l2 < graphY) y2_l2 = graphY;

            DrawLineEx((Vector2){x1, y1_l1}, (Vector2){x2, y2_l1}, 2.0f, MAGENTA);
            DrawLineEx((Vector2){x1, y1_l2}, (Vector2){x2, y2_l2}, 2.0f, CYAN);
        }
        DrawText("Eigenvalue History", graphX, graphY - 20, 10, DARKGRAY);

        EndDrawing();
    }
    CloseWindow();
    return 0;
}
