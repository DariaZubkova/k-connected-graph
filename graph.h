#pragma once
#include "matrix.h"
#include <map>
#include <vector>
#include <list>

const int INF = 1000000000;

class Graph {
public:
	Graph() {}
	~Graph() {}

	int DFS(std::list <std::pair <int, std::pair <int, int> > >& block, std::vector <std::vector <int> >& lay, int& n, int& s, int& t, int cost);
	void build(std::vector <std::vector <int> >& lay, int& n, int& s, int& t);
	int Dinitz(int& s, int& t);
	bool algorithmGalil();
	bool connect(int s);


	Graph(int numNode_, int k_);
	bool algorithmEven();
	void extraMatrix(int numI, int numJ);
	Matrix get_Matrix();
	Matrix get_ExtraMatrix();
	bool checkMinGraph();
	bool checkContractionMinmality();
	bool perebor_vertex();
	bool perebor_vertex_fun(int &num, std::vector<int> &index);
	bool perebor_vertex_check();
	int get_k();
	int get_numNode();
	int get_numEdge();
	int get_numLine();
	int get_numColumn();
	int get_num();
	void set_Matrix(Matrix m);
	void set_ExtraMatrix(Matrix m);
	void set_k(int k_);
	void set_numNode(int nN);
	void set_numEdge(int nE);
	void set_numLine(int nL);
	void set_numColumn(int nC);
	void set_num(int n);
	bool bfs(Matrix f, std::vector<int>& d, int s, int t);
	int dfs(int u, int minFlow, Matrix& f, std::vector<int> ptr, std::vector<int> d, int t);
	int dinic(int s, int t);
	void set_minimality(bool res);
	void set_minContraction(bool res);
	void addEdge(int J, int pos_s);
private:
	Matrix matrix;
	Matrix extra_matrix;
	int numNode;
	int numEdge;
	int k;
	int numLine;
	int numColumn;
	int num;
	bool minimality;
	bool minContraction;
};
