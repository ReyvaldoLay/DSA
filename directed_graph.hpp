#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include<iostream>
#include<string>
#include<vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <stack>
#include <limits>
#include <utility> 
#include <algorithm>
// include more libraries here if you need to

using namespace std; // the standard namespace are here just in case.

/*
	the vertex class
*/
template <typename T>
class vertex {

public:
	int id;
	T weight;

	vertex(int x, T y) : id(x), weight(y) {}

	// add more functions here if you need to
};

/*
	the graph class
*/
template <typename T>
class directed_graph {

private:

	//You will need to add some data members here
	//to actually represent the graph internally,
	//and keep track of whatever you need to.
	unordered_map<int, T> vertex_weights;
	unordered_map<int, unordered_map<int, T>> adj_list;

public:

	directed_graph(); //A constructor for directed_graph. The graph should start empty.
	~directed_graph(); //A destructor. Depending on how you do things, this may not be necessary.

	bool contains(const int&) const; //Returns true if the graph contains the given vertex_id, false otherwise.
	bool adjacent(const int&, const int&); //Returns true if the first vertex is adjacent to the second, false otherwise.

	void add_vertex(const vertex<T>&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const int&, const int&, const T&); //Adds a weighted edge from the first vertex to the second.

	void remove_vertex(const int&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const int&, const int&); //Removes the edge between the two vertices, if it exists.

	size_t in_degree(const int&) ; //Returns number of edges coming in to a vertex.
	size_t out_degree(const int&) ; //Returns the number of edges leaving a vertex.
	size_t degree(const int&) ; //Returns the degree of the vertex (both in edges and out edges).

	size_t num_vertices() const; //Returns the total number of vertices in the graph.
	size_t num_edges() const; //Returns the total number of edges in the graph.

	vector<vertex<T>> get_vertices(); //Returns a vector containing all the vertices.
	vector<vertex<T>> get_neighbours(const int&); //Returns a vector containing all the vertices reachable from the given vertex. The vertex is not considered a neighbour of itself.
	vector<vertex<T>> get_second_order_neighbours(const int&); // Returns a vector containing all the second_order_neighbours (i.e., neighbours of neighbours) of the given vertex.
															  // A vector cannot be considered a second_order_neighbour of itself.
	bool reachable(const int&, const int&) ; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
	bool contain_cycles() ; // Return true if the graph contains cycles (there is a path from any vertices directly/indirectly to itself), false otherwise.

	vector<vertex<T>> depth_first(const int&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
	vector<vertex<T>> breadth_first(const int&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

	directed_graph<T> out_tree(const int&); //Returns a tree starting at the given vertex using the out-edges. This means every vertex in the tree is reachable from the root.

	vector<vertex<T>> pre_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of a pre-order traversal of the tree starting at the given vertex.
	vector<vertex<T>> in_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of an in-order traversal of the tree starting at the given vertex.
	vector<vertex<T>> post_order_traversal(const int&, directed_graph<T>&); // returns the vertices in ther visitig order of a post-order traversal of the tree starting at the given vertex.

	vector<vertex<T>> significance_sorting(); // Return a vector containing a sorted list of the vertices in descending order of their significance.
};

// Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
// Although these are just the same names copied from above, you may find a few more clues in the full method headers.
// Note also that C++ is sensitive to the order you declare and define things in - you have to have it available before you use it.

template <typename T>
directed_graph<T>::directed_graph() {}

template <typename T>
directed_graph<T>::~directed_graph() {}

template <typename T>
bool directed_graph<T>::contains(const int& u_id) const{ 
	if(vertex_weights.find(u_id)!=vertex_weights.end()){ // checks vertex_weights if it has a vertex with that id
		return true;
	}
	return false;
}

template <typename T>
bool directed_graph<T>::adjacent(const int& u_id, const int& v_id){
	if(contains(u_id) && contains(v_id)){ //checks if both are in the graph
		if(adj_list[u_id].find(v_id)!=adj_list[u_id].end()){  //checks the adj_list of u if it has a vertex of v_id
			return true;
		}
	}
	return false;
}

template <typename T>
void directed_graph<T>::add_vertex(const vertex<T>& u) {
	if(!contains(u.id)){ //checks if it is already in the graph | only pass when it isnt in the graph 
		vertex_weights.insert({u.id, u.weight}); // add to vertex_weights 
		adj_list[u.id]=unordered_map<int, T>(); // create a adj_list for the vertex
	}
}

template <typename T>
void directed_graph<T>::add_edge(const int& u_id, const int& v_id, const T& weight) {
	if(contains(u_id) && contains(v_id)){ //checks that both are in graph | only pass if both are in graph
		if(adj_list[u_id].find(v_id)==adj_list[u_id].end()){ // checks if v exists in adj_list of u | pass when doesnt
			adj_list[u_id].insert({v_id, weight}); // add v to adj_list of u with weight
		}
	}
}

template <typename T>
void directed_graph<T>::remove_vertex(const int& u_id) {
	if(contains(u_id)){ //checks if it exists in graph | pass when it does
		vertex_weights.erase(u_id); // looks through vertex_weights and erase the vertex
		adj_list.erase(u_id); // look through adj_list and erase the list for the vertex
		for (auto& x: adj_list){ // iterate through adj_list
			x.second.erase(u_id); // checks the second element and erase it if it is the vertex
		}
	}
}

template <typename T>
void directed_graph<T>::remove_edge(const int& u_id, const int& v_id) {
	if(adjacent(u_id, v_id)){ // checks if there is an edge from u to v
		adj_list[u_id].erase(v_id); // go to the adj_list of u and remove v
	}
}

template <typename T>
size_t directed_graph<T>::in_degree(const int& u_id)  {
	size_t in_num = 0; // initialize the thing to return | 0 so it never returns null
	for(auto& x: adj_list){ // iterate through every member of adj_list
		if(adjacent(x.first, u_id)){ // checks is u is in the adj_list of every member | pass if u is in the current adj_list
			in_num++; // increase the counter by 1
		}
	}
	return in_num;
}

template <typename T>
size_t directed_graph<T>::out_degree(const int& u_id){
	size_t out_num = 0; // initialize the thing to return | 0 so it never returns null
	for(auto& x: adj_list[u_id]){ // iterate through every member of adj_list[u]
		out_num++; // increase the counter by 1
	}
	return out_num;
}

template <typename T>
size_t directed_graph<T>::degree(const int& u_id) {
	size_t deg_num= out_degree(u_id) + in_degree(u_id); // calls other degree function to not write extra code | even if both dont point to or pointed to it will return 0 not null
	return deg_num;
 }

template <typename T>
size_t directed_graph<T>::num_vertices() const {
	size_t count = 0; //initialize thing to return | 0 so it never returns null 
	for (auto& x : vertex_weights){ //iterate trough every vertex
		count ++; // increase the counter by 1
	}
	return count;
}

template <typename T>
size_t directed_graph<T>::num_edges() const {
	size_t count = 0; //initialize thing to return | 0 so it never returns null
	for (auto& x: adj_list){ // iterate through all adj_list
		count += x.second.size(); // adds the number of vertex the current list is pointing at to the counter
	}
	return count;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_vertices() {
	vector<vertex<T>> v; // initialize thing to return so it can be added to
	for(auto x: vertex_weights){ // iterate through all vertex
		v.push_back(vertex<T>(x.first, x.second)); // add the vertex to the local vector
	}
	return v;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_neighbours(const int& u_id) {
	vector<vertex<T>> v; // initialize thing to return so it can be added to
	if(contains(u_id)){ // checks if u is in the graph
		for (auto x: adj_list[u_id]){ // iterate through all members of adj_list[u] | adj_list[u] only contains things that it is pointing to
			v.insert(v.begin(),vertex<T>(x.first, vertex_weights[x.first])); // add the vertex to the local vector
		}
	}
	return v;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_second_order_neighbours(const int& u_id) {
	vector<vertex<T>> v; // // initialize thing to return so it can be added to
	if(contains(u_id)){ // checks if u exists in the graph
		for(auto x: get_neighbours(u_id)){  // iterate through all vertex that u is pointing to
			for(auto y: get_neighbours(x.id)){ // iterate through every member of the neighbours
				if(y.id != u_id){ // is it isnt the same as what was inputed | pass when different
					bool is_duplicate = false; // initialize boolean value to be used later
					for(auto z : v){ // iterate through every member of the local vector
    					if (z.id == y.id){ // if the y is already in the local vector
       						is_duplicate = true; // set to true so it doesnt get added again
        					break; // get out of the for loop
    					 }
					}
					if (!is_duplicate){ // if itsnt a duplicate
    					v.push_back(vertex(y.id, vertex_weights[y.id])); //add it to local vector
					}
				}
			}
		}
	}
	return v;
}

template <typename T>
bool directed_graph<T>::reachable(const int& u_id, const int& v_id){
	bool *visited = new bool[num_vertices()]; // create an array of boolean to check if it is visited or not
	for (unsigned i = 0; i < num_vertices(); i++){ // mark all as false or not visited
		visited[i] = false;
	}
	
	queue<int> unprocessed; // initialize the queue
	unprocessed.push(u_id); // put the first element in
	
	vector<vertex<T>> ordered; // initialize the vector to add vertexes to
	
	while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.front(); // get the front element
		unprocessed.pop(); // remove the front element
		ordered.push_back(vertex<T>(n, vertex_weights[n])); // add the element with the same number as the front element to local vector

		if(get_neighbours(n).size() > 0){ // checks if it has any neighbours
			for(auto x : get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					visited[x.id] = true; // set it as visited 
					unprocessed.push(x.id); // add the queue so it can be processed last
				}
			}
		}
	}
	ordered.erase(ordered.begin()); //remove the first element which is always the same as u_id
	for(auto x : ordered){ // iterate through ordered
		if(x.id == v_id){ // checks if the id of the current vertex is the same as v_id | pass if same
			return true; 
		}
	}
	return false;
}

template <typename T>
bool directed_graph<T>::contain_cycles()  {
	for(auto x : get_vertices()){ // iterate through all vertices
		if(reachable(x.id, x.id)){ // checks if it is reachable to itself | pass if it is
			return true;
		}
	}
	return false;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::depth_first(const int& u_id) {
	bool *visited = new bool[num_vertices()]; // create an array of boolean to check if it is visited or not
	for (unsigned i = 0; i < num_vertices(); i++){ // mark all as false or not visited
		visited[i] = false;
	}
	
	stack<int> unprocessed; // initialize the stack
	unprocessed.push(u_id); // put the first element in
	
	vector<vertex<T>> ordered; // initialize the vector to add vertexes to
	
	while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.top();  // get the top element
		visited[n] = true; //set visited as true so it doesnt process it again
		unprocessed.pop(); // remove the top element
		ordered.push_back(vertex<T>(n, vertex_weights[n])); // add the element with the same number as the front element to local vector

		if(get_neighbours(n).size() > 0){ // checks if it has any neighbours
			for(auto x : get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					visited[x.id] = true; // set as visited
					unprocessed.push(x.id); // add to top of stack so it is processed next
				}
			}
		}
	}
	return ordered;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::breadth_first(const int& u_id) {
	bool *visited = new bool[num_vertices()]; // create an array of boolean to check if it is visited or not
	for (unsigned i = 0; i < num_vertices(); i++){  // mark all as false or not visited
		visited[i] = false;
	}
	
	queue<int> unprocessed; // initialize the queue
	unprocessed.push(u_id); // put the first element in
	
	vector<vertex<T>> ordered; // initialize the vector to add vertexes to
	
	while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.front(); // get the front element
		visited[n] = true; //set visited as true so it doesnt process it again
		unprocessed.pop(); // remove the front element
		ordered.push_back(vertex<T>(n, vertex_weights[n])); // add the element with the same number as the front element to local vector

		if(get_neighbours(n).size() > 0){ // checks if it has any neighbours
			for(auto x : get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					visited[x.id] = true; // set as visited
					unprocessed.push(x.id); // add to the back of queue so it is processed last
				}
			}
		}
	}
	return ordered;
}

template <typename T>
directed_graph<T> directed_graph<T>::out_tree(const int& u_id) {
	directed_graph<T> g; // initialize thing to return and so stuff can be added to it
    vector<bool> visited(num_vertices(), false);
    queue<pair<int, int>> q;

    if(u_id >= 0){
        q.push({u_id, u_id});
        while (!q.empty()){
            pair<int, int> n = q.front();
            q.pop();
            if (!visited[n.first]){
                visited[n.first] = true;
                if(n.first != n.second){
                    g.add_vertex(vertex<T>(n.first, vertex_weights[n.first]));
                    g.add_vertex(vertex<T>(n.second, vertex_weights[n.second]));
                    g.add_edge(n.second, n.first, 0);
                }
                for (auto i : get_vertices()){
                    if (adj_list[n.first][i.id]){
                        q.push({i.id, n.first});
                    }
                }
            }
        }
    }
    return g;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::pre_order_traversal(const int& u_id, directed_graph<T>& mst) {
	bool *visited = new bool[mst.num_vertices()]; // create an array of boolean to check if it is visited or not
	for (unsigned i = 0; i < mst.num_vertices(); i++){ // mark all as false or not visited
		visited[i] = false;
	}
	
	stack<int> unprocessed; // initialize the stack
	unprocessed.push(u_id); // put the first element in
	
	vector<vertex<T>> preOrder; // initialize the vector to add vertexes to
	
	while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.top();  // get the top element
		visited[n] = true; //set visited as true so it doesnt process it again
		unprocessed.pop(); // remove the top element
		preOrder.push_back(vertex<T>(n, mst.vertex_weights[n])); // add the element with the same number as the front element to local vector

		if(mst.get_neighbours(n).size() > 0){ // checks if it has any neighbours
			for(auto x : mst.get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					visited[x.id] = true; // set as visited
					unprocessed.push(x.id); // add to top of stack so it is processed next
				}
			}
		}
	}
	return preOrder;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::in_order_traversal(const int& u_id, directed_graph<T>& mst) {
	vector<vertex<T>> inOrder;
	// bool *visited = new bool[num_vertices()]; // create an array of boolean to check if it is visited or not
	// for (unsigned i = 0; i < num_vertices(); i++){  // mark all as false or not visited
	// 	visited[i] = false;
	// }
	// stack<int> unprocessed;
	// unprocessed.push(u_id);
	// vector<vertex<T>> nb;
	// while(!unprocessed.empty()){
	// 	int counter = 0;
	// 	int n = unprocessed.top();
	// 	visited[n] = true;

	// 		nb = mst.get_neighbours(n);
	// 		for(auto x : nb){
	// 			if(!visited[x.id]){
	// 				unprocessed.push(x.id);
	// 				counter++;
	// 				break;
	// 			}
	// 		}
	// 		if(counter == 0 && nb.size() == 0){
	// 			unprocessed.pop();
	// 			inOrder.push_back(vertex<T>(n, mst.vertex_weights[n]));

	// 		}
	// 		if(counter == 0 && nb.size() > 0){
	// 			unprocessed.pop();
    // 			inOrder.push_back(vertex(n, mst.vertex_weights[n]));
	// 		}
	// }

// 	stack<int> s; 
//     int n = u_id;
//   vector<vertex<T>> nb;
// 	while (mst.contains(n) || s.empty() == false) {
// 		while (mst.contains(n)) {
// 			s.push(n); 
// 			nb = mst.get_neighbours(n);
//             n = nb[0].id; 
//         }
// 		n = s.top(); 
//         s.pop();
// 		inOrder.push_back(vertex<T>(n, vertex_weights[n]));
// 		nb = mst.get_neighbours(n);
//         n = nb[1].id; 
// 	}

// 	vector<vertex<T>> child = mst.get_neighbours(u_id);
// 	vertex<T> right = child[0];
// 	vertex<T> left = child[1];
// 	if(u_id ==0){
// 		return inOrder;
// 	}
// 	in_order_traversal(left.id, mst);
// 	inOrder.push_back(vertex<T>(u_id, vertex_weights[u_id]));
// 	in_order_traversal(right.id, mst);

	return inOrder;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::post_order_traversal(const int& u_id, directed_graph<T>& mst) {
	vector<vertex<T>> postOrder;
	bool *visited = new bool[num_vertices()]; // create an array of boolean to check if it is visited or not
	for (unsigned i = 0; i < num_vertices(); i++){  // mark all as false or not visited
		visited[i] = false;
	}
	
	stack<int> unprocessed; // initialize the queue
	unprocessed.push(u_id); // put the first element in
	
	while (!unprocessed.empty()){ // loops when unprocessed has stuff in it
		int n = unprocessed.top(); // get the top
		visited[n] = true; // set visited to true so it doesnt process it again
		int counter = 0; // 

			for(auto x : mst.get_neighbours(n)){ // iterate though all neighbours of n
				if(!visited[x.id]){ // check if it has been visited | pass if it has not
					unprocessed.push(x.id); // add to the top of stack so it is processed next
					counter++; // add to counter so if counter = 0 it doesnt have any neighbours that hasnt been visited
				}
			}
		if(mst.get_neighbours(n).size() == 0 || counter == 0){ //  if it doesnt have any neighbours at all(a leaf) || if all neighbours are already visited/processed
			unprocessed.pop(); // remove from last so it is not processed again
			postOrder.push_back(vertex<T>(n, vertex_weights[n])); // add a vertex with the same id to the final vector
		}
	}
	return postOrder;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::significance_sorting() { // using bubble sort
	vector<vertex<T>> s = get_vertices();
		bool swapped; // can terminate early if no swap occurs
	do{
		swapped = false;
		for(int i = 0; i < s.size()-1; i++){
			if(s[i].weight > s[i+1].weight){ // compare its weight with the vertex and swap them if the next is smaller.
				vertex<T> temp = s[i];
				s[i] = s[i+1];
				s[i+1] = temp;
				swapped = true; 
			}
		} // always move the largest to the right each iteration
	}while(swapped); 
	return s;
}

#endif
