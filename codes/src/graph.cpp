#include "src/graph.h"

Graph::Graph(int num_vertices, VertexInfo * vertices_info)
{
	num_vertices_ = num_vertices;
    arcs_ = new Matrix<GArc*>(num_vertices,num_vertices, nullptr);
	num_arcs_ = 0;

	adj_lists_out_ = new std::vector<std::list<int>>(num_vertices);
    adj_lists_in_ = new std::vector<std::list<int>>(num_vertices);

	vertices_info_ = vertices_info;
	indexes_ = new Matrix<int>(num_vertices,num_vertices,-1);
}

Graph::~Graph()
{
    if(vertices_info_)
    {
        delete [] vertices_info_;
        vertices_info_ = nullptr;
    }

    if(indexes_)
    {
        delete indexes_;
        indexes_ = nullptr;
    }

    if(adj_lists_out_)
    {
        for(int v1 = 0; v1 < num_vertices_; ++v1)
        {
            for(std::list<int>::iterator it = ((*adj_lists_out_)[v1]).begin(); it != ((*adj_lists_out_)[v1]).end(); ++it)
            {
                int v2 = *it;
                delete (*arcs_)[v1][v2];
                (*arcs_)[v1][v2] = nullptr;
            }
        }
        delete adj_lists_out_;
        adj_lists_out_ = nullptr;
    }

    if(adj_lists_in_)
    {
        delete adj_lists_in_;
        adj_lists_in_ = nullptr;
    }

    if(arcs_)
    {
        delete arcs_;
        arcs_ = nullptr;
    }
}

GArc** Graph::operator[](int i) const
{
	return (*arcs_)[i];
}

std::ostream& operator<< (std::ostream &out, const Graph &g)
{
	GArc * arc = nullptr;
    for(int i=0; i < g.num_vertices(); ++i)
	    out << i << ") " <<(g.vertices_info())[i].coordinates_.first << " " <<(g.vertices_info())[i].coordinates_.second << " " <<(g.vertices_info())[i].profit_  << " " << (g.vertices_info())[i].decay_ratio_
            << " " << (g.vertices_info())[i].nominal_service_time_ << " " << (g.vertices_info())[i].dev_service_time_ << std::endl;

    for(int i=0; i < g.num_vertices(); ++i)
        for(int j=0; j < g.num_vertices(); ++j)
			if(arc = (g[i][j]))
                out << "(" << i << "," << j << ") " <<  *arc << std::endl;
                
    return out;
}

void Graph::set_num_vertices(int num)
{
	num_vertices_ = num;
}

void Graph::set_indexes(Matrix<int> * indexes)
{
    if(indexes_)
    {
        delete indexes_;
        indexes_ = nullptr;
    }
    arcs_map_.clear();
    indexes_ = indexes;
}

void Graph::SetArcIndex(int v1, int v2, int index)
{
    (*indexes_)[v1][v2] = index;
    arcs_map_[index] = std::pair<int,int>(v1,v2);
}

void Graph::AddArc(int i, int j, double dist)
{
    GArc * curr_arc = new GArc(dist);
    (*arcs_)[i][j] = curr_arc;
    ((*adj_lists_out_)[i]).push_back(j);
    ((*adj_lists_in_)[j]).push_back(i);
    (*indexes_)[i][j] = num_arcs_;
    arcs_map_[num_arcs_] = std::pair<int,int>(i,j);
    ++num_arcs_;
}

void Graph::DeleteArc(int i, int j, bool remove_from_indexes_list)
{
    if(remove_from_indexes_list)
    {
        ((*adj_lists_out_)[i]).remove(j);
        ((*adj_lists_in_)[j]).remove(i);
    }
    delete ((*arcs_)[i][j]);
    ((*arcs_)[i][j]) = nullptr;
    arcs_map_[(*indexes_)[i][j]] = std::pair<int,int>(-1,-1);
    (*indexes_)[i][j] = -1;
    --num_arcs_;
}

void Graph::AddEdge(int i, int j, double dist)
{
    AddArc(i,j,dist);
    AddArc(j,i,dist);
}

int Graph::pos(int i, int j) const
{
    return (*indexes_)[i][j];
}

void Graph::WriteGraph(std::string file_name) const
{
	GArc * curr_arc = NULL;
	std::fstream file;
	file.open(file_name.c_str(),std::fstream::out);

	file << num_vertices_ << std::endl;

    for(int i=0; i < num_vertices_; ++i)
    {
        for(int j=0; j < num_vertices_; ++j)
        {
            curr_arc = (*arcs_)[i][j];
            if(curr_arc)
                file << i << " " << j << " " << curr_arc->distance() << std::endl;
        }
    }
}
