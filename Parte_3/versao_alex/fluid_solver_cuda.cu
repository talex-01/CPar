#include <cuda_runtime.h>
#include <stdio.h>

__global__ void red_black_step(int M, int N, int O, float *x, const float *x0, float a, float c, bool is_red) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i > M || j > N || k > O) return;

    // Determine if the current cell is red or black
    bool is_current_red = (i + j + k) % 2 == 0;

    if (is_current_red != is_red) return;

    int index = i + (M + 2) * (j + (N + 2) * k);
    int xm = index - 1;
    int xp = index + 1;
    int ym = index - (M + 2);
    int yp = index + (M + 2);
    int zm = index - (M + 2) * (N + 2);
    int zp = index + (M + 2) * (N + 2);

    x[index] = (x0[index] + a * (x[xm] + x[xp] + x[ym] + x[yp] + x[zm] + x[zp])) / c;
}

void lin_solve_cuda(int M, int N, int O, int b, float *x, float *x0, float a, float c) {
    size_t size = (M + 2) * (N + 2) * (O + 2) * sizeof(float);
    float *d_x, *d_x0;

    // Allocate memory on the GPU
    cudaMalloc(&d_x, size);
    cudaMalloc(&d_x0, size);

    // Copy data to GPU
    cudaMemcpy(d_x, x, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_x0, x0, size, cudaMemcpyHostToDevice);

    // Define CUDA grid and block dimensions
    dim3 blockDim(8, 8, 8);
    dim3 gridDim((M + 7) / 8, (N + 7) / 8, (O + 7) / 8);

    const int max_iters = 20;
    for (int iter = 0; iter < max_iters; ++iter) {
        // Red pass
        red_black_step<<<gridDim, blockDim>>>(M, N, O, d_x, d_x0, a, c, true);
        cudaDeviceSynchronize();

        // Black pass
        red_black_step<<<gridDim, blockDim>>>(M, N, O, d_x, d_x0, a, c, false);
        cudaDeviceSynchronize();
    }

    // Copy results back to CPU
    cudaMemcpy(x, d_x, size, cudaMemcpyDeviceToHost);

    // Free GPU memory
    cudaFree(d_x);
    cudaFree(d_x0);
}
