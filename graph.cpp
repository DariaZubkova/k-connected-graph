#include"graph.h"
#include <mutex>
#include <fstream>
std::mutex mtx_dinic;
#include <queue>

Graph::Graph(int numNode_, int k_) {
	numNode = numNode_;
	k = k_;
	numLine = numNode;
	numColumn = numNode;
	Matrix matrix_(numLine, numColumn);
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			if (i == j) {
				matrix_.setElem(i, j, 0);
			}
			else {
				matrix_.setElem(i, j, 1);
			}
		}
	}
	matrix = matrix_;
	num = 0;
	minimality = false;
	minContraction = false;
}

bool Graph::bfs(Matrix f, std::vector<int>& d, int s, int t) {
	int numLine = extra_matrix.getNumLine();
	int numColumn = extra_matrix.getNumColumn();
	if (numLine == 0 || numColumn == 0) {
		return false;
	}
	std::vector<int> Q;
	int u = 0;
	for (int i = 0; i < numLine; i++) {
		d[i] = -1;
	}
	d[s] = 0;
	Q.push_back(s);
	while (!Q.empty()) {
		u = Q[0];
		Q.erase(Q.begin());
		for (int v = 0; v < numLine; v++) {
			int value_f = f.getElem(u, v);
			int value_extra = extra_matrix.getElem(u, v);
			if (d[v] == -1 && value_f < value_extra) { //&& f.getElem(u, v) < extra_matrix.getElem(u, v)
				Q.push_back(v);
				d[v] = d[u] + 1;
			}
		}
	}
	Q.clear();
	return d[t] != -1;
}

int min(int a, int b) {
	if (a > b)
		return b;
	else
		return a;
}

Graph extraGraph(Graph graph, int J) {
	int numLine = graph.get_numLine(), numColumn = graph.get_numColumn();
	if (numLine == 0 || numColumn == 0 || graph.get_numNode() <= 0) {
		Graph newGraph;
		return newGraph;
	}
	int newNumLine = numLine + 1, newNumColumn = numColumn + 1;
	Matrix newMatrix(newNumLine, newNumColumn);
	Matrix matrix = graph.get_Matrix();
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			int value = matrix.getElem(i, j);
			newMatrix.setElem(i, j, value);
		}
	}
	int i = newNumLine - 1;
	int j = 0;
	while (j < J) {
		newMatrix.setElem(i, j, 1);
		newMatrix.setElem(j, i, 1);
		j++;
	}
	Graph newGraph;
	newGraph.set_k(graph.get_k());
	newGraph.set_num(graph.get_num());
	newGraph.set_numNode(graph.get_numNode() + 1);
	newGraph.set_Matrix(newMatrix);
	newGraph.set_numLine(newNumLine);
	newGraph.set_numColumn(newNumColumn);
	newGraph.set_numEdge(newMatrix.currentNumEdge());
	return newGraph;
}

void Graph::addEdge(int J, int pos_s) {
	matrix.setElem(J, pos_s, 1);
	matrix.setElem(pos_s, J, 1);
}

int Graph::dfs(int u, int minFlow, Matrix& f, std::vector<int> ptr, std::vector<int> d, std::vector<int> &vertex, int t) {
	int numLine = extra_matrix.getNumLine();
	int numColumn = extra_matrix.getNumColumn();
	if (numLine == 0 || numColumn == 0) {
		return 0;
	}
	if (u == t || minFlow == 0) {
		return minFlow;
	}
	for (int v = ptr[u]; v < numLine; v++) {
		int value_f = f.getElem(u, v);
		int value_extra = extra_matrix.getElem(u, v);
		if (d[v] == d[u] + 1) { //&& value_extra != 0  && value_extra > 0
			if (std::find(vertex.begin(), vertex.end(), v) == vertex.end()) {
				int delta = dfs(v, min(minFlow, value_extra - value_f), f, ptr, d, vertex, t);
				if (delta > 0) {
					vertex.push_back(u);
					f.addVal(u, v, delta);
					f.difVal(v, u, delta);
					extra_matrix.difVal(u, v, delta);
					extra_matrix.addVal(v, u, delta);
					return delta;
				}
			}
		}
	}
	return 0;
}

int Graph::dinic(int s, int t) {
	int numLine = extra_matrix.getNumLine();
	int numColumn = extra_matrix.getNumColumn();
	int maxFlow = 0;
	int flow = 0;
	if (numLine == 0 || numColumn == 0) {
		return 0;
	}
	Matrix f(numLine, numColumn);
	for (int h = 0; h < numLine; h++) {
		for (int l = 0; l < numColumn; l++) {
			f.setElem(h, l, 0);
		}
	}
	std::vector<int> ptr(numLine, 0);
	std::vector<int> vertex;
	std::vector<int> d(numLine, -1);
	while (bfs(f, d, s, t)) {
		for (int i = 0; i < numLine; i++)
			ptr[i] = 0;
		flow = dfs(s, INF, f, ptr, d, vertex, t);
		if (flow == 0) {
			return maxFlow;
		}
		while (flow != 0) {
			maxFlow += flow;
			if (maxFlow >= k) {
				return maxFlow;
			}
			flow = dfs(s, INF, f, ptr, d, vertex, t);
		}
	}
	return maxFlow;
}

bool Graph::algorithmEven() {
	int maxFlow = -1;
	if (numLine <= 0 || numColumn <= 0) {
		return false;
	}
	for (int j = 1; j < k; j++) {
		for (int i = 0; i < j; i++) {
			(*this).extraMatrix(i, j);
			maxFlow = (*this).dinic(i, j);
			if (maxFlow < k) {
				return false;
			}
		}
	}
	int pos_s = matrix.getNumColumn();
	Graph newGraph = extraGraph((*this), k);
	for (int j = k; j < matrix.getNumColumn(); j++) {
		newGraph.extraMatrix(j, pos_s);
		maxFlow = newGraph.dinic(j, pos_s);
		if (maxFlow < k) {
			return false;
		}
		newGraph.addEdge(j, pos_s);
	}
	return true;
}

bool Graph::algorithmGalil() {
	int maxFlow = -1;
	if (numLine <= 0 || numColumn <= 0) {
		return false;
	}
	int kk = 1, MIN = numLine - 1;
	while (true) {
		std::vector<int> maxFlowVec;
		for (int i = 0; i < kk - 1; i++) {
			(*this).extraMatrix(i, kk);
			maxFlow = (*this).dinic(i, kk);
			maxFlowVec.push_back(maxFlow);
		}
		int minFlow = 1000000;
		for (size_t i = 0; i < maxFlowVec.size(); i++) {
			if (minFlow > maxFlowVec[i]) {
				minFlow = maxFlowVec[i];
			}
		}
		MIN = std::min(minFlow, MIN);
		if (!maxFlowVec.empty()) {
			maxFlowVec.clear();
		}
		if (MIN < kk) {
			std::cout << "MIN connected = " << MIN << std::endl;
			if (MIN == k)
				return true;
			else
				return false;
		}
		int pos_s = matrix.getNumColumn();
		Graph newGraph = extraGraph((*this), k);
		for (int j = kk; j < matrix.getNumColumn(); j++) {
			newGraph.extraMatrix(j, pos_s);
			maxFlow = newGraph.dinic(j, pos_s);
			if (maxFlow < kk) {
				std::cout << "G connected = " << kk - 1 << std::endl;
				if (kk - 1 == k)
					return true;
				else
					return false;
			}
			newGraph.addEdge(j, pos_s);
		}
		kk++;
	}
	return true;
}

void Graph::extraMatrix(int numI, int numJ) {
	int extraVertex = 0;
	std::map<int, int> n;
	if (numLine <= 0 || numColumn <= 0) {
		Matrix newMatrix(0, 0);
		extra_matrix = newMatrix;
		return ;
	}
	std::vector<int> dVertex(numLine, 0);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine; i++) {
		if (i != numI && i != numJ) {
			std::pair<int, int> extraNode;
			extraNode.first = i;
			extraNode.second = numLine + extraVertex;
			n.insert(extraNode);
			extraNode.first = numLine + extraVertex;
			extraNode.second = i;
			n.insert(extraNode);
			extraVertex++;
		}
	}
	int newNumLine = numLine + extraVertex;
	int newNumColumn = numColumn + extraVertex;
	Matrix newMatrix(newNumLine, newNumColumn);
	for (int i = 0; i < newNumLine; i++) {
		for (int j = 0; j < newNumColumn; j++) {
			if (i != numI && i != numJ) {
				if (i < numLine) {
					auto prevJ = n.find(i);
					if (prevJ != n.end()) {
						if (j == prevJ->second) {
							newMatrix.setElem(i, j, dVertex[i]); //dVertex[i]
						}
						else {
							newMatrix.setElem(i, j, 0);
						}
					}
				}
				else {
					if (j < numColumn) {
						auto prevI = n.find(i);
						if (prevI != n.end()) {
							int value = matrix.getElem(prevI->second, j);
							newMatrix.setElem(i, j, value);
							//newArr[newNumColumn * i + j] = arr[numColumn * prevI->second + j];//map
						}
					}
					else {
						newMatrix.setElem(i, j, 0);
					}
				}
			}
			else {
				if (j < numColumn) {
					int value = matrix.getElem(i, j);
					newMatrix.setElem(i, j, value);
					//newArr[newNumColumn * i + j] = arr[numColumn * i + j];
				}
				else {
					newMatrix.setElem(i, j, 0);
				}
			}
		}
	}
	n.clear();
	extra_matrix = newMatrix;
}

Matrix Graph::get_Matrix() {
	Matrix m(numLine, numColumn);
	m = matrix;
	return m;
}

Matrix Graph::get_ExtraMatrix() {
	Matrix m;
	m = extra_matrix;
	return m;
}

void Graph::set_minimality(bool res) {
	minimality = res;
}

bool Graph::checkMinGraph() {
	int maxFlow = 0;
	std::vector<int> dVertex(numLine, 0);
	if (numLine <= 0 || numColumn <= 0) {
		return false;
	}
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
					(*this).extraMatrix(i, j);
					maxFlow = (*this).dinic(i, j);
					if (maxFlow >= k) {
						return false;
					}
					matrix.setElem(j, i, 1);
				}
				matrix.setElem(i, j, 1);
			}
		}
	}
	return true;
}

bool Graph::checkContractionMinmality() {
	bool res = false;
	int maxFlow = 0;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1) {
				Graph new_graph(numNode - 2, k - 1); //k - 1
				Matrix new_matrix(numNode - 2, numNode - 2);
				int pos_i = 0, pos_j = 0;
				for (int k = 0; k < numLine; k++) {
					if (k != i && k != j) {
						for (int l = 0; l < numColumn; l++) {
							if (l != i && l != j) {
								new_matrix.setElem(pos_i, pos_j, matrix.getElem(k, l));
								pos_j++;
							}
						}
						pos_j = 0;
						pos_i++;
					}
				}
				new_graph.set_Matrix(new_matrix);

				//res = new_graph.algorithmGalil();
				res = new_graph.algorithmEven();
				if (res == true) {
					return false;
				}
			}
		}
	}
	return true;
}

void Graph::set_minContraction(bool res) {
	minContraction = res;
}

int Graph::get_k() {
	if (k <= 0) {
		return 0;
	}
	int nK = k;
	return nK;
}

int Graph::get_numNode() {
	if (numNode <= 0) {
		return 0;
	}
	int nN = numNode;
	return nN;
}

int Graph::get_numEdge() {
	int nE = matrix.currentNumEdge();
	return nE;
}

int Graph::get_numLine() {
	int nL = matrix.getNumLine();
	return nL;
}

int Graph::get_numColumn() {
	int nC = matrix.getNumColumn();
	return nC;
}

int Graph::get_num() {
	if (num <= 0) {
		return 0;
	}
	int n = num;
	return n;
}

void Graph::set_Matrix(Matrix m) {
	matrix = m;
}

void Graph::set_ExtraMatrix(Matrix m) {
	extra_matrix = m;
}

void Graph::set_k(int k_) {
	k = k_;
}

void Graph::set_numNode(int nN) {
	numNode = nN;
}

void Graph::set_numEdge(int nE) {
	numEdge = nE;
}

void Graph::set_numLine(int nL) {
	numLine = nL;
}

void Graph::set_numColumn(int nC) {
	numColumn = nC;
}

void Graph::set_num(int n) {
	num = n;
}