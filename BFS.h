#include <algorithm>
#include <queue>
#include <set>
#include <cmath>
#include <climits>

#include "parallel.h"

using namespace parlay;
using namespace std;
#define THRESHOLD 10000
#define FLT_SIZES_TH 100
#define FLT_TH 1000

template<typename T>
void pscan(T* A, size_t n) {
  if(n==0) return;
  if(n<=THRESHOLD) {
    for(size_t i=1;i<n;i++) A[i] += A[i-1];
    return;
  }
  //From homework 2 prb 5
  // Step 1 Split the array into k chunks, each of size about ⌈n/k⌉ except for probably the last one
  size_t k = ceil(sqrt(n));
  size_t k_size = n/k;
  if(n%k != 0) k_size++;

  //Step 2  sums of k chucks are stored in an array S [1..k]
  T* S = new T[k];
  S[0] = 0;

  //Compute sum of chunks in parallel
  parallel_for(1, k, [&](size_t i) { 
    S[i] = 0;
    for(size_t j = (i-1)*k_size;j<i*k_size;j++) {
      if(j>=n) break;
      S[i] += A[j];
    }
  });
  //Step 3 prefix sum of S[1..k] sequentially
  for(size_t i=1;i<k;i++) S[i] += S[i-1];

  //Step 4 Create k parallel task for k chunks and 
  //compute the prefix sum of a chunk based on input array 
  //and S calculated in previous step
   parallel_for(0, k, [&](size_t i) { 
    for(size_t j = 0;j<k_size;j++) {
      size_t ind = j+k_size*i;
      if(ind < n) {
        if(j==0) A[ind] += S[i];
        else A[ind] += A[ind-1];
      }
    }
  });
  delete[] S;
}

//Packing
uint32_t* pack(uint64_t* offsets, uint32_t* edges, uint32_t* flag, size_t &n, size_t node) {
  //No. of neighbors for current node
  size_t n_ngh = offsets[node+1] - offsets[node];
  uint32_t* B = new uint32_t[n_ngh];
  pscan(flag,n_ngh);
  n = flag[n_ngh-1];
  //output array to contains the node neighbors
  if (flag[0] == 1) B[0] = edges[offsets[node]];
  // Coarsening
  if(n_ngh < THRESHOLD) {
    for(size_t i = 1;i < n_ngh;i++) {
      if (flag[i] != flag[i-1])
                B[flag[i] - 1] = edges[i+offsets[node]];
    }
  } else {
    parallel_for( 1,  n_ngh, [&](size_t i) { 
      if (flag[i] != flag[i-1])
                B[flag[i] - 1] = edges[i+offsets[node]];
  });
  }

  return B;
}

uint32_t* flatten(uint32_t** A, uint32_t* S, size_t &n){
    uint32_t* offset = new uint32_t[n];
    parallel_for(0, n, [&](size_t i) { offset[i] = S[i];});
    // std::cout<<"\n-----------Inside Flatten-----\n";
    //  for(size_t i=0;i<n;i++) {
    //   std::cout<<"S[i] ::"<<S[i];
    //   for(uint32_t j=0;j<S[i];j++) std::cout<<A[i][j]<<" ";
    //   std::cout<<endl;
    // }
    pscan(offset,n);
    size_t new_n = offset[n-1];
    uint32_t* B = new uint32_t[new_n];
    if(n < FLT_TH) {
      for(size_t i=0;i<n;i++){
        //Since the scan is inclusive
        size_t off = i!=0?offset[i-1]:0;
        if(S[i] < FLT_SIZES_TH) {
          for(uint32_t ngh = 0;ngh<S[i];ngh++) B[off + ngh] = A[i][ngh];
        } else {
          parallel_for(0, S[i],[&](uint32_t ngh) { B[off + ngh] = A[i][ngh];});
        }
        delete[] A[i];
      }
    } else {
      parallel_for(0,n,[&](size_t i){
        size_t off = i!=0?offset[i-1]:0;
        if(S[i] < FLT_SIZES_TH) {
          for(uint32_t ngh = 0;ngh<S[i];ngh++) B[off + ngh] = A[i][ngh];
        } else {
          parallel_for(0, S[i],[&](uint32_t ngh) { B[off + ngh] = A[i][ngh];});
        }
        delete[] A[i];
      });

    }
    delete[] offset;
    delete[] A;
    delete[] S;
    n = new_n;
    return B;
}
uint32_t* new_frontier_sparse(uint64_t* offsets, uint32_t* edges, uint32_t* dist, uint32_t* curr_frontier, size_t &f_n) {
  uint32_t** out_nghs = new uint32_t*[f_n];
  uint32_t* S = new uint32_t[f_n];

  
  parallel_for(0,f_n,[&](size_t i) {

    size_t node = curr_frontier[i];
    uint32_t off = offsets[node];
    size_t n_ngh = offsets[node+1] - offsets[node];
        size_t nghs_size = 0;

    uint32_t* out_flag = new uint32_t[n_ngh];
    if(n_ngh<FLT_TH) {
      for(size_t j=0;j<n_ngh;j++) out_flag[j] =0;
      for(uint32_t j = offsets[node]; j<offsets[node +1];j++) {
        uint32_t ngh = edges[j];
        if(dist[ngh] == INT_MAX && __sync_bool_compare_and_swap(&dist[ngh], INT_MAX, (dist[node] +1)))
          out_flag[j-off] = 1;
      }
    } else {
        parallel_for(0,n_ngh,[&](size_t j) { out_flag[j] =0;});
        parallel_for(offsets[node],offsets[node +1],[&](uint32_t j){
          uint32_t ngh = edges[j];
          if(dist[ngh] == INT_MAX && __sync_bool_compare_and_swap(&dist[ngh], INT_MAX, (dist[node] +1)))
            out_flag[j-off] = 1;
      });
    }
    // std::cout<<"\n------Out flag-----\n";
    // for(size_t j =0;j<n_ngh;j++) std::cout<<out_flag[j]<<std::endl;
    //     std::cout<<"\n------Inside Distance-----\n";

    //     for(size_t j =0;j<6;j++) std::cout<<dist[j]<<std::endl;

    out_nghs[i] = pack(offsets,edges,out_flag,nghs_size,node);
    S[i] = nghs_size;
    delete[] out_flag;
  });
  return flatten(out_nghs,S,f_n);
}
void BFS(uint64_t* offsets, uint32_t* edges, uint32_t* dist, size_t n, size_t m,
         uint32_t s) {
  parallel_for(0,n,[&](size_t i) {dist[i] = INT_MAX;});
  dist[s] = 0;
  size_t f_n = 1;
  uint32_t* frontier = new uint32_t[f_n];
  frontier[0] = s;
  //std::cout<<"\n--------Frontier-----------\n";
  while(f_n > 0) {
    frontier = new_frontier_sparse(offsets, edges,dist,frontier,f_n);
  }
  delete[] frontier;
}
