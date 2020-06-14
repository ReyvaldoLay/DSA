#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

using namespace std;

/*
 * Computes the shortest distance from u to v in graph g.
 * The shortest path corresponds to a sequence of vertices starting from u and ends at v,
 * which has the smallest total weight of edges among all possible paths from u to v.
 */
template <typename T>
vector<vertex<T>> shortest_path(directed_graph<T> g, int u_id, int v_id) {
  vector<vertex<T>> neighbours;
  vector<vertex<T>> vertices;
  set<pair<int,int>> setdest; 

  vector<T> dist(g.num_vertices()+1, 9999999); //create distance array to store shortest distance
  vector<int> prev(g.num_vertices()+1, 0); //create preview array to store previous shortest vertex
  
  setdest.insert(make_pair(0, u_id)); //insert source and set dest to 0
  
  dist[u_id] = 0; 
  
  while(!setdest.empty()) { 
    pair<int, int> temp = *(setdest.begin()); //make the top vertex in set the one we are working with
    setdest.erase(setdest.begin()); //remove it from set so that it isnt processed again
    
    int u = temp.second; //set temp vertex to vertex id
    
    neighbours = g.get_neighbours(u); //get neighbours of the current vertex
  
    for(auto x: neighbours) { 
      int v = x.id; //Neighbour id
      int weight = g.edge_weight(u, v); //get current-neighbour edge cost
      
      if(dist[v] > dist[u] + weight) { //if current distance is smaller than total then go into code 
        if(dist[v] != 9999999) {
          setdest.erase(setdest.find(make_pair(dist[v], v))); //update dist[neighbour] with new distance
        }
        dist[v] = dist[u] + weight; //set new distance
        prev[v] = u; //set previous of shortest neighbour to current
        setdest.insert(make_pair(dist[v], v)); //add to set
      }
    }
  }
  //trace back from v_id to u_id
  int currentVertex = v_id; //Set current to target
  
  while(prev[currentVertex] != 0) { //while previous is not 0
    vertices.push_back(g.get_vertex(currentVertex)); //add the vertex to final list 
    currentVertex = prev[currentVertex]; //set current to prev
  } 
  
  vertices.push_back(g.get_vertex(u_id)); //add u_id since it is the source/ start
  reverse(vertices.begin(), vertices.end()); //reverse the vector for final vector
  
  return vertices;
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */
template <typename T>
vector<vector<vertex<T>>> strongly_connected_components(directed_graph<T> g) {
  vector<vector<vertex<T>>> scc;
  int visited[g.num_vertices()+1] = {0};
  int low[g.num_vertices()+1];
  int visited_count = 0;
  stack<vertex<T>> st;
  set<int> se; // visited vertex id
  
  for(auto x : g.get_vertices()) { //loop through all vertices
    if(!visited[x.id]) { // if unvisitied then call scc_util
      SCC_Util(x, visited, low, visited_count, st, se, scc, g);
    }
  }
  return scc;
}

template <typename T>
void SCC_Util(vertex<T>& x, int visited[], int low[], int& visited_count, stack<vertex<T>>& st, set<int>& se, vector<vector<vertex<T>>>& scc, directed_graph<T>& g) {
  low[x.id] = visited[x.id] = ++visited_count; // add to visited and increment visited count
  st.push(x);
  se.insert(x.id);
  
  for(auto y : g.get_neighbours(x.id)) { // for each neighbour of current vertex
    // check if the subtree rooted with 'y' has a way to get to 'x'
    if(!visited[y.id]) { // if unvisited
      SCC_Util(y, visited, low, visited_count, st, se, scc, g); // recursive
      low[x.id] = min(low[x.id], low[y.id]);
    }
    else if(se.find(y.id) != se.end()) { // update low value of 'x' if 'y' is still in the stack
      low[x.id] = min(low[x.id], visited[y.id]);
    }
  }
  
  if(visited[x.id] == low[x.id]) { // check if current vertex is root vertex
    vector<vertex<T>> vertices; // create vector of to keep vertices
    
    
    while(st.top().id != x.id) {//add all vertices above current vertex to vector
      vertices.push_back(st.top());
      se.erase(st.top().id); // remove from vertex id
      st.pop(); // remove from stack
    }
    vertices.push_back(st.top()); // add x to list
    se.erase(st.top().id); // remove from vertex id
    st.pop(); // remove from stack
    scc.push_back(vertices); // add vertices to scc
  }
}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 * You will be given a DAG as the argument.
 */
template <typename T>
vector<vertex<T>> topological_sort(directed_graph<T> g) {
  vector<vertex<T>> sorted;
  vector<bool> visited(g.num_vertices(), false);
  queue<int> unprocessed; // initialize the queue
  for(auto x : g.get_vertices()){
    if(g.in_degree(x.id) == 0){ //get all starting vertices
      unprocessed.push(x.id); // add to queue
    }
  }
  while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.front(); // get the front element
		visited[n] = true; //set visited as true so it doesnt process it again
		unprocessed.pop(); // remove the front element
		sorted.push_back(g.get_vertex(n)); // add the element with the same number as the front element to local vector

		if(g.get_neighbours(n).size() > 0){ // checks if it has any neighbours
			for(auto x : g.get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					visited[x.id] = true; // set as visited
					unprocessed.push(x.id); // add to the back of queue so it is processed last
				}
			}
		}
	}
  return sorted;
}

/*
 * Computes the lowest cost-per-person for delivery over the graph.
 * u is the source vertex, which send deliveries to all other vertices.
 * vertices denote cities; vertex weights denote cities' population;
 * edge weights denote the fixed delivery cost between cities, which is irrelevant to 
 * the amount of goods being delivered. 
 */
template <typename T>
T low_cost_delivery(directed_graph<T> g, int u_id) {
  size_t cost = 0;
  size_t dividend;
  size_t divisor;
  for(auto x : g.get_vertices()){ // get all weight of vertices
    divisor = divisor + x.weight; // add all weight of vertices to divisor
  }
  divisor = divisor-g.get_vertex(u_id).weight; // remove the initial/ starting points weight from divisor

  vector<bool> visited(g.num_vertices(), false); // create a vector with the initial value of false
  queue<int> unprocessed; // initialize the queue
  unprocessed.push(u_id); // add to queue

  while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.front(); // get the front element
		visited[n] = true; //set visited as true so it doesnt process it again
		unprocessed.pop(); // remove the front element

		if(g.get_neighbours(n).size() > 0){ // checks if it has any neighbours
			for(auto x : g.get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					visited[x.id] = true; // set as visited
					unprocessed.push(x.id); // add to the back of queue so it is processed last
          dividend = dividend + g.edge_weight(n,x.id); // add the weight of edge so it can be processed later
				}
			}
		}
	}
  cost = dividend / divisor;
  return cost;
}

/*
i also added 2 functions to directed_graph.hpp in the second assignment
template <typename T>
size_t directed_graph<T>::edge_weight(const int& u_id, const int& v_id) {
	size_t count = 0;
	for(auto x : adj_list[u_id]){
		if(x.first == v_id){
			return(x.second);
		}
	}
	return count;
}

template <typename T>
vertex<T> directed_graph<T>::get_vertex(const int& u_id) {
	vertex<T> v = vertex<T>(0,0);
	for(auto x: vertex_weights){
		if(x.first == u_id){
			v = vertex<T>(x.first,x.second);
		}
	}
	return v;
}
*/
