#pragma once
#include"matrix.h"
const int INF = 1000000000;

class Graph {
public:
	Graph() {}
	~Graph() {}
	Graph(int numNode_, int k_);
	bool algorithmEven();
	void extraMatrix(int numI, int numJ);
	Matrix get_Matrix();
	Matrix get_ExtraMatrix();
	bool checkMinGraph();
	bool checkContractionMinmality();
	int get_k();
	int get_numNode();
	int get_numEdge();
	int get_numLine();
	int get_numColumn();
	int get_num();
	void set_Matrix(Matrix& m);
	void set_ExtraMatrix(Matrix& m);
	void set_k(int k_);
	void set_numNode(int nN);
	void set_numEdge(int nE);
	void set_numLine(int nL);
	void set_numColumn(int nC);
	void set_num(int n);
	bool bfs(Matrix& f, std::vector<int>& d, int s, int t);
	int dfs(int u, int minFlow, Matrix& f, std::vector<int>& ptr, std::vector<int>& d, int s, int t);
	int dinic(Matrix f, std::vector<int> ptr, std::vector<int> d, int s, int t);
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