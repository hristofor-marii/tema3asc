#include <iostream>
#include <stdlib.h>
#include "helper.h"
#include <fstream>
#include <cuda.h>
#include <cuda_runtime_api.h>

vector<Input> take_in(const char* fileIn)
{
    string geon1;
    float lat1;
    float lon1;
    int pop1;
    vector<Input> in;
    ifstream ifs(fileIn);

    while(ifs >> geon1 >> lat1 >> lon1 >> pop1) {

        Input aux;
        aux.geon = geon1;
        aux.lat = lat1;
        aux.lon = lon1;
        aux.pop = pop1;

        in.push_back(aux);
    }

    ifs.close();

    return in;
}

__device__ float geoDistances(float lat1, float lon1, float lat2, float lon2)
{
	float phi1 = (90.f - lat1) * DEGREE_TO_RADIANS;
    	float phi2 = (90.f - lat2) * DEGREE_TO_RADIANS;
    	float theta1 = lon1 * DEGREE_TO_RADIANS;
    	float theta2 = lon2 * DEGREE_TO_RADIANS;
	float cs = sin(phi1) * sin(phi2) * cos(theta1 - theta2) + cos(phi1) * cos(phi2);
    
	if (cs > 1) {
       		cs = 1;
    	} else if (cs < -1) {
       		cs = -1;
    	}

    	return 6371.f * acos(cs);
}
__global__ void calc_pop(const float *lat, const float *lon, const int *pop, int *acc_pop, const size_t km_range, const size_t n)
{
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if(i < n && j < n) {
		if(i < j) {

			float distance = geoDistances(lat[i], lon[i], lat[j], lon[j]);
			if(distance <= km_range)
			{
				atomicAdd(&acc_pop[i], pop[j]);
				atomicAdd(&acc_pop[j], pop[i]);
			}
                }
                else if(i == j)
                {
                        atomicAdd(&acc_pop[i], pop[j]);
                }
        }

}

int main(int argc, char** argv) {
    DIE( argc == 1,
         "./accpop <kmrange1> <file1in> <file1out> ...");
    DIE( (argc - 1) % 3 != 0,
         "./accpop <kmrange1> <file1in> <file1out> ...");

    for(int argcID = 1; argcID < argc; argcID += 3) {
        float kmRange = atof(argv[argcID]);
        vector<Input> in;

	in = take_in(argv[argcID + 1]);	       
	ofstream ofs(argv[argcID + 2]);
	float *lat = 0;
	float *lon = 0;
	int *pop = 0;
	int *acc_pop = 0;
	int num_elem = in.size();

	cudaMallocManaged(&lat, num_elem * sizeof(float));
	cudaMallocManaged(&lon, num_elem * sizeof(float));
	cudaMallocManaged(&pop, num_elem * sizeof(int));
	cudaMallocManaged(&acc_pop, num_elem * sizeof(int));
	
	if (lat == 0 || lon == 0 || pop == 0 || acc_pop == 0) {
        	cout << "[HOST] Couldn't allocate memory\n";
        	return 1;
    	}
	
	for(int i = 0; i < num_elem; i++)
       	{
		lat[i] = in[i].lat;
		lon[i] = in[i].lon;
		pop[i] = in[i].pop;
		acc_pop[i] = 0;
	}
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	dim3 block(32,32);
	dim3 grid((num_elem + 31)/32,  (num_elem + 31)/32);

	cudaEventRecord(start);
	calc_pop<<<grid, block>>>(lat, lon, pop, acc_pop, kmRange, num_elem);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	/*
       	float ms = 0;
	cudaEventElapsedTime(&ms, start, stop);
	float seconds = ms / pow((float) 10, 3);

	cout<< "time = " << seconds  << endl;
	*/
	for(int i = 0; i < num_elem; i++)
		ofs << acc_pop[i] << endl;

	cudaFree(lat);
	cudaFree(lon);
	cudaFree(pop);
	cudaFree(acc_pop);
    }
}
