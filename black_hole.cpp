#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <fstream>
#include <sstream>
#ifdef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;
using Clock = std::chrono::high_resolution_clock;

// VARS
double lastPrintTime = 0.0;
int frameCount = 0;
double c = 299792458.0; // speed of light in m/s
double G = 6.67430e-11; // gravitational constant in m^3 kg^-1 s^-2
struct Ray;
bool Gravity = true;

struct Camera {
    // Center the camera orbit on the black hole at (0,0,0)
    vec3 target = vec3(0.0f, 0.0f, 0.0f); // Always look at the black hole center
    float radius = 6.34194e10f; // Distance from the black hole (in meters)
    float minRadius = 1e10f, maxRadius = 1e12f; // Zoom limits

    float azimuth = 0.0f; // Horizontal angle
    float elevation = M_PI / 2.0f; // Vertical angle

    float orbitSpeed = 0.01f; // Speed of orbiting around the black hole
    float panSpeed = 0.01f; // Speed of panning the camera
    float zoomSpeed = 1.1f; // Speed of zooming in/out

    bool dragging = false;
    bool panning = false;
    bool moving = false;
    double lastX = 0.0, lastY = 0.0;

    //Calculate the camera's position in Cartesian coordinates
    vec3 position() const {
        float clampedElevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        // Orbiting around (0,0,0)
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update(){
        // Always keep target at the black hole center
        target = vec3(0.0f, 0.0f, 0.0f);
        if(dragging | panning){
            moving = true;
        } else {
            moving = false;
        }
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);

        if(dragging && panning) {
            // Pan: Shift + Left or Middle button
            // Disable panning to keep camera centered on black hole
        }
        else if(dragging && !panning) {
            // Orbit: Left button
            azimuth += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f); // Prevent flipping
        }

        lastX = x;
        lastY = y;
        update();
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if(button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if(action == GLFW_PRESS) {
                dragging = true;
                // Disable panning to keep camera centered on black hole
                panning = false;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if(action == GLFW_RELEASE) {
                dragging = false;
                panning = false;
            }
        }
        if(button == GLFW_MOUSE_BUTTON_RIGHT) {
            if(action == GLFW_PRESS) {
                Gravity = true;
            } else if(action == GLFW_RELEASE) {
                Gravity = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        radius = glm::clamp(radius, minRadius, maxRadius);
        update();
    }
    void processKey(int key, int scancode, int action, int mods) {
        if(action == GLFW_PRESS && key == GLFW_KEY_G) {
            Gravity = !Gravity;
            cout << "[INFO] Gravity turned " << (Gravity ? "ON" : "OFF") << endl;
        }
    }
};
Camera camera;

struct BlackHole {
    vec3 position;
    double mass;
    double radius;
    double r_s;

    BlackHole(vec3 pos, float m) : position(pos), mass(m) {r_s = 2.0 * G * mass / (c*c);}
    bool Intercept(float px, float py, float pz) const {
        double dx = double(px) - double(position.x);
        double dy = double(py) - double(position.y);
        double dz = double(pz) - double(position.z);
        double dist2 = dx*dx + dy*dy + dz*dz;
        return dist2 < r_s * r_s;
    }
};
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36); // Mass of Sagittarius A black hole
struct ObjectData {
    vec4 posRadius; // xyz = position, w = radius
    vec4 color;    // rgb = color, a = unused
    float mass;
    vec3 velocity = vec3(0.0f, 0.0f, 0.0f); // Initial velocity
};
vector<ObjectData> objects = {
    {vec4(4e11f, 0.0f, 0.0f, 1e10f), vec4(1,1,0,1), 1.98891e30 },
    {vec4(0.0f, 0.0f, 4e11f, 4e10f), vec4(1,0,0,1), 1.98891e30 },
    {vec4(0.0f, 0.0f, 0.0f, SagA.r_s), vec4(0,0,0,1), static_cast<float>(SagA.mass) },
    // {vec4(6e10f, 0.0f, 0.0f, 5e10f), vec4(0,1,0,1) }
};