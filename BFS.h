#include <algorithm>
#include <queue>
#include <set>

#include "parallel.h"

using namespace parlay;
using namespace std;

void BFS(uint64_t* offsets, uint32_t* edges, uint32_t* dist, size_t n, size_t m,
         uint32_t s) {
  dist[s] = 0;
  queue<uint32_t> q;
  q.push(s);
  while (!q.empty()) {
    uint32_t u = q.front();
    q.pop();
    for (size_t i = offsets[u]; i < offsets[u + 1]; i++) {
      uint32_t v = edges[i];
      if (dist[v] > dist[u] + 1) {
        dist[v] = dist[u] + 1;
        q.push(v);
      }
    }
  }
}
