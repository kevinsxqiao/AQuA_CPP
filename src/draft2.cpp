#include "data/data.h"

std::vector<std::vector<Point4D>> bwconncomp4D(std::vector<std::vector<std::vector<std::vector<bool>>>>& BW) {
    std::vector<std::vector<Point4D>> components;
    std::vector<std::vector<std::vector<std::vector<bool>>>> visited(BW.size(), std::vector<std::vector<std::vector<bool>>>(BW[0].size(), std::vector<std::vector<bool>>(BW[0][0].size(), std::vector<bool>(BW[0][0][0].size(), false))));

    int dx[3] = {-1, 0, 1};
    int dy[3] = {-1, 0, 1};
    int dz[3] = {-1, 0, 1};
    int dt[3] = {-1, 0, 1};

    for (int t = 0; t < BW.size(); ++t) {
        for (int z = 0; z < BW[0].size(); ++z) {
            for (int x = 0; x < BW[0][0].size(); ++x) {
                for (int y = 0; y < BW[0][0][0].size(); ++y) {
                    if (BW[t][z][x][y] && !visited[t][z][x][y]) {
                        std::vector<Point4D> component;
                        std::queue<Point4D> queue;
                        queue.push(Point4D(x, y, z, t));
                        visited[t][z][x][y] = true;

                        while (!queue.empty()) {
                            Point4D current = queue.front();
                            queue.pop();
                            component.push_back(current);

                            for (int dti = 0; dti < 3; ++dti) {
                                for (int dzi = 0; dzi < 3; ++dzi) {
                                    for (int dxi = 0; dxi < 3; ++dxi) {
                                        for (int dyi = 0; dyi < 3; ++dyi) {
                                            int nx = current.x + dx[dxi];
                                            int ny = current.y + dy[dyi];
                                            int nz = current.z + dz[dzi];
                                            int nt = current.t + dt[dti];

                                            if (nx >= 0 && nx < BW[0][0].size() && ny >= 0 && ny < BW[0][0][0].size() && nz >= 0 && nz < BW[0].size() && nt >= 0 && nt < BW.size() && BW[nt][nz][nx][ny] && !visited[nt][nz][nx][ny]) {
                                                queue.push(Point4D(nx, ny, nz, nt));
                                                visited[nt][nz][nx][ny] = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        components.push_back(component);
                    }
                }
            }
        }
    }

    return components;
}
