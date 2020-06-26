#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <cassert>
#include <chrono>
#include <math.h>

using namespace std;

float lb = 0;
float ub = 2;
int nx = 11;
int ny = 11;
int nt = 100;
int nit = 50;
int c = 1;
float dx = ub / float(nx - 1);
float dy = ub / float(ny - 1);


float rho = 1;
float nu = 0.1;
float dt = 0.001;

void linspace(vector<float> &x, int lb, int ub, int num)
{
    for (int k = 0; k < num; k++)
    {
        x[k] = k * (ub - lb) / float(num - 1) + lb;
    }
}

void linspace(float* x, int lb, int ub, int num)
{
    for (int k = 0; k < num; k++)
    {
        x[k] = k * (ub - lb) / float(num - 1) + lb;
    }
}

void meshgrid(vector<float> x, vector<float> y,
              vector<vector<float>> &X, vector<vector<float>> &Y)
{
    int dimx = x.size();
    int dimy = y.size();
    assert(dimy == X.size() && dimy == Y.size());

    for(int k=0; k< dimy; k++)
        X[k] = x;
    for(int k=0; k< dimy; k++)
        Y[k] = vector<float>(dimx, y[k]); 

}

void build_up_b(vector<vector<float>> &b, float rho, float dt,
                vector<vector<float>> u, vector<vector<float>> v,
                float dx, float dy)
{
    /*

def build_up_b(b, rho, dt, u, v, dx, dy):
    return b
*/
    int xlim = b.size();
    int ylim = b[0].size();

    for (int idx = 0; idx < xlim - 2; idx++)
        for (int idy = 0; idy < ylim - 2; idy++)
        {
            b[idx + 1][idy + 1] = (rho * (1 / dt *
                                              ((u[idx + 1][idy + 2] - u[idx + 1][idy]) /
                                                   (2 * dx) +
                                               (v[idx + 2][idy + 1] - v[idx][idy + 1]) / (2 * dy)) -
                                          pow(((u[idx + 1][idy + 2] - u[idx + 1][idy]) / (2 * dx)), 2) -
                                          2 * ((u[idx + 2][idy + 1] - u[idx][idy + 1]) / (2 * dy) *
                                               (v[idx + 1][idy + 2] - v[idx + 1][idy]) / (2 * dx)) -
                                          pow((v[idx + 2][idy + 1] - v[idx][idy + 1]) / (2 * dy), 2)));
        }
}

void pressure_poisson(vector<vector<float>> &p, float dx, float dy, vector<vector<float>> b)
{
    /*
def pressure_poisson(p, dx, dy, b):

*/
    vector<vector<float>> pn = p;
    int xlim = p.size();
    int ylim = p[0].size();

    for (int q = 0; q < nit; q++)
    {
        for (int idx = 0; idx < xlim - 2; idx++)
            for (int idy = 0; idy < ylim - 2; idy++)
            {
                pn = p;
                p[idx + 1][idy + 1] = (((pn[idx + 1][idy + 2] + pn[idx + 1][idy]) * pow(dy, 2) +
                                        (pn[idx + 2][idy + 1] + pn[idx][idy + 1]) * pow(dx, 2)) /
                                           (2 * (pow(dx, 2) + pow(dy, 2))) -
                                       pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
                                           b[idx + 1][idy + 1]);
            }
        p[0] = p[1];                          // dp/dy = 0 at y = 0
        p[xlim - 1] = vector<float>(ylim, 0); // p = 0 at y = 2

        for (int idx = 0; idx < xlim; idx++)
        {
            p[idx][ylim - 1] = p[idx][ylim - 2]; // dp/dx = 0 at x = 2
            p[idx][0] = p[idx][1];               // dp/dx = 0 at x = 0
        }
    }
}

void cavity_flow(int nt, vector<vector<float>> &u, vector<vector<float>> &v,
                 float dt, float dx, float dy, vector<vector<float>> &p, float rho, float nu)
{
    /*
    def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
*/
    int xlim = u.size();
    int ylim = u[0].size();
    vector<vector<float>> un = u;
    vector<vector<float>> vn = v;
    vector<vector<float>> b(ny, vector<float>(nx));

    for (int n = 0; n < nt; n++)
    {
        un = u;
        vn = v;
        build_up_b(b, rho, dt, u, v, dx, dy);
        pressure_poisson(p, dx, dy, b);

        for (int idx = 0; idx < xlim - 2; idx++)
            for (int idy = 0; idy < ylim - 2; idy++)
            {
                u[idx + 1][idy + 1] = (un[idx + 1][idy + 1] -
                                       un[idx + 1][idy + 1] * dt / dx *
                                           (un[idx + 1][idy + 1] - un[idx + 1][idy]) -
                                       vn[idx + 1][idy + 1] * dt / dy *
                                           (un[idx + 1][idy + 1] - un[idx][idy + 1]) -
                                       dt / (2 * rho * dx) * (p[idx + 1][idy + 2] - p[idx + 1][idy]) +
                                       nu * (dt / pow(dx, 2) *
                                                 (un[idx + 1][idy + 2] - 2 * un[idx + 1][idy + 1] + un[idx + 1][idy]) +
                                             dt / pow(dy, 2) *
                                                 (un[idx + 2][idy + 1] - 2 * un[idx + 1][idy + 1] + un[idx][idy + 1])));

                v[idx + 1][idy + 1] = (vn[idx + 1][idy + 1] -
                                       un[idx + 1][idy + 1] * dt / dx *
                                           (vn[idx + 1][idy + 1] - vn[idx + 1][idy]) -
                                       vn[idx + 1][idy + 1] * dt / dy *
                                           (vn[idx + 1][idy + 1] - vn[idx][idy + 1]) -
                                       dt / (2 * rho * dy) * (p[idx + 2][idy + 1] - p[idx][idy + 1]) +
                                       nu * (dt / pow(dx, 2) *
                                                 (vn[idx + 1][idy + 2] - 2 * vn[idx + 1][idy + 1] + vn[idx + 1][idy]) +
                                             dt / pow(dy, 2) *
                                                 (vn[idx + 2][idy + 1] - 2 * vn[idx + 1][idy + 1] + vn[idx][idy + 1])));
            }
        u[0] = vector<float>(ylim, 0);
        u[xlim - 1] = vector<float>(ylim, 10);
        v[0] = vector<float>(ylim, 0);
        v[xlim - 1] = vector<float>(ylim, 0);
        for (int idx = 0; idx < xlim; idx++)
        {
            u[idx][0] = 0;
            v[idx][0] = 0;
            u[idx][ylim - 1] = 0;
            v[idx][ylim - 1] = 0;
        }
    }
}

void npprint(vector<vector<float>> u)
{
    cout << "[" << endl;
    for (int i = 0; i < u.size(); i++)
    {
        cout << "[";
        for (int k = 0; k < u[0].size(); k++)
            cout << u[i][k] << ", ";
        cout << "]," << endl;
    }
    cout << "]" << endl;
    cout << "x-------------------------------x" << endl;
}

int main()
{
    // vector<float> x(nx);
    float x[nx];
    vector<float> y(ny);
    vector<vector<float>> X(ny);
    vector<vector<float>> Y(ny);
    linspace(x, lb, ub, nx);
    npprint(x);
    // linspace(y, lb, ub, ny);
    // meshgrid(x, y, X, Y);
    // npprint(X);
    // npprint(Y);
    // cout << x.size() << endl;

    vector<vector<float>> u(ny, vector<float>(nx));
    vector<vector<float>> v(ny, vector<float>(nx));
    vector<vector<float>> b(ny, vector<float>(nx));
    vector<vector<float>> p(ny, vector<float>(nx));

    auto start = std::chrono::high_resolution_clock::now();

    cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu);
    auto finish = std::chrono::high_resolution_clock::now();
    // npprint(u);
    // npprint(v);
    // npprint(p);

    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

}