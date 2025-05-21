#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include "src/arc.h"
#include "src/matrix.hpp"

class Graph final{
public:
	struct VertexInfo{
		std::pair<double,double> coordinates_;
		int profit_;
		double decay_ratio_;
		double nominal_service_time_;
		double dev_service_time_;
	};

	explicit Graph() = default;
	explicit Graph(int,VertexInfo*);
	virtual ~Graph();
	int num_vertices() const {return num_vertices_;}
	int num_arcs(void) const {return num_arcs_;}
	std::list<int>& AdjVerticesOut(int vertex) const{
		return (*adj_lists_out_)[vertex];
	}

	std::list<int>& AdjVerticesIn(int vertex) const{
		return (*adj_lists_in_)[vertex];
	}

	VertexInfo * vertices_info() const {return vertices_info_;}
	void set_num_vertices(int);
	void set_indexes(Matrix<int> * indexes);
	void SetArcIndex(int v1, int v2, int index);
	int pos(int i, int j) const;
	void AddArc(int, int, double);
	void DeleteArc(int, int, bool remove_from_indexes_lists = true);
	void AddEdge(int, int, double);
	GArc** operator[](int) const;
	void WriteGraph(std::string) const;
	friend std::ostream& operator<< (std::ostream &out, const Graph &graph);
private:
	Matrix<GArc*> * arcs_ = nullptr;
	VertexInfo * vertices_info_ = nullptr;
	int num_vertices_ = 0;
	int num_arcs_ = 0;
	std::vector<std::list<int>> * adj_lists_out_ = nullptr;
	std::vector<std::list<int>> * adj_lists_in_ = nullptr;
	Matrix<int> * indexes_ = nullptr;
	std::unordered_map<int,std::pair<int,int>> arcs_map_;
};
