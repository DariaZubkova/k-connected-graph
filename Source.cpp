#include "matrix.h"
#include"graph.h"
#include <map>
#include <queue>
#include <utility>
#include <string>
#include <ctime>
#include <fstream>

bool compareMatrix(std::vector<Graph>& arr_graph, Matrix& matrix) {
	int nL = matrix.getNumLine(), nC = matrix.getNumColumn();
	bool res = true;
	//matrix.PrintMatrix();
	//std::cout << "COMPARE:" << std::endl;
	for (auto v : arr_graph) {
		res = true;
		//v.PrintMatrix();
		//std::cout << std::endl;
		Matrix m = v.get_Matrix();
		for (int i = 0; i < nL; i++) {
			for (int j = 0; j < nC; j++) {
				if (m.getElem(i, j) != matrix.getElem(i, j)) {
					res = false;
				}
			}
		}
		if (res) {
			return res;
		}
	}
	//std::cout << std::endl;
	return res;
}

void fun(Graph& graph, int* num, int k, std::vector<Graph>* arr_graph) { //перебор графов
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	int maxFlow;
	std::vector<int> dVertex;
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex.push_back(s);
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				graph.set_Matrix(matrix);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				res = graph.algorithmEven();
				if (res && !compareMatrix((*arr_graph), matrix)) {
					//std::cout << *num + 1 << std::endl;
					(*arr_graph).push_back(graph);
				}
				(*num) += 1;
				fun(graph, num, k, arr_graph);
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}

void check_k_connected(std::vector<Graph> arr) {
	int numVector = 0;
	if (!arr.empty()) {
		int k = arr[0].get_k();
		for (auto v : arr) {
			std::cout << "Matrix number: " << numVector << std::endl;
			Matrix m = v.get_Matrix();
			m.PrintMatrix();
			int numEdge = m.currentNumEdge();
			std::cout << "numEdge " << numEdge << std::endl;
			std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
			std::cout << std::endl;
			numVector++;
		}
	}
	std::cout << "END!" << std::endl;
}

void check() {
	int numNode = 100;
	int numEdge = 0;
	int k = 7;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);
	unsigned int start_time = clock();
	bool res = false;
	res = graph.algorithmEven();
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	res = graph.checkMinGraph();
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - end_time;
	std::cout << "Time work = " << search_time2 / CLOCKS_PER_SEC << std::endl;
	res = graph.checkContractionMinmality();
	unsigned int end_time3 = clock();
	unsigned int search_time3 = end_time3 - end_time2;
	std::cout << "Time work = " << search_time3 / CLOCKS_PER_SEC << std::endl;
}

void printToFile(std::vector<Graph> graphs, std::string nameFile) {
	std::fstream f;
	f.open(nameFile, std::ios::out);
	if (f)
	{
		std::string str;

		for (auto f_graph : graphs)
		{
			int numLine = f_graph.get_numLine();
			int numColumn = f_graph.get_numColumn();
			Matrix m = f_graph.get_Matrix();
			for (int i = 0; i < numLine; i++) {
				str = "T ";
				str += std::to_string(i);
				f << str << std::endl;
				str = "";
			}
			for (int i = 0; i < numLine; i++) {
				for (int j = i+1; j < numColumn; j++) {
					if (m.getElem(i, j) == 1) {
						str = "E ";
						str += std::to_string(i);
						str += " ";
						str += std::to_string(j);
						str += " ";
						str += std::to_string(m.getElem(i, j));
						f << str << std::endl;
						str = "";
					}
				}
			}
		}
		f.close();
	}
}

int main(void) {
	int numNode = 6;
	int numEdge = 0;
	int k = 4;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);
	unsigned int start_time = clock();
	std::cout << "Start algorithm " << std::endl;
	std::cout << std::endl;
	bool res;
	int num = 0;
	std::vector<Graph> arr_graph;
	
	/*
	Matrix matrixCheck(numLine + 5, numColumn + 5);
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			if (i == j) {
				matrixCheck.setElem(i, j, 0);
				//arr[numColumn * i + j] = 0;
			}
			else {
				matrixCheck.setElem(i, j, 1);
				//arr[numColumn * i + j] = 1;
			}
		}
	}
	matrixCheck.setElem(0, 4, 1);
	matrixCheck.setElem(0, 5, 1);
	matrixCheck.setElem(0, 6, 1);
	matrixCheck.setElem(0, 7, 1);
	matrixCheck.setElem(0, 8, 1);
	matrixCheck.setElem(1, 4, 1);
	matrixCheck.setElem(1, 5, 1);
	matrixCheck.setElem(1, 6, 1);
	matrixCheck.setElem(1, 7, 1);
	matrixCheck.setElem(1, 8, 1);
	matrixCheck.setElem(2, 4, 1);
	matrixCheck.setElem(2, 5, 1);
	matrixCheck.setElem(2, 6, 1);
	matrixCheck.setElem(2, 7, 1);
	matrixCheck.setElem(2, 8, 1);
	matrixCheck.setElem(3, 4, 1);
	matrixCheck.setElem(3, 5, 1);
	matrixCheck.setElem(3, 6, 1);
	matrixCheck.setElem(3, 7, 1);
	matrixCheck.setElem(3, 8, 1);
	matrixCheck.setElem(4, 0, 1);
	matrixCheck.setElem(4, 1, 1);
	matrixCheck.setElem(4, 2, 1);
	matrixCheck.setElem(4, 3, 1);
	matrixCheck.setElem(5, 0, 1);
	matrixCheck.setElem(5, 1, 1);
	matrixCheck.setElem(5, 2, 1);
	matrixCheck.setElem(5, 3, 1);
	matrixCheck.setElem(6, 0, 1);
	matrixCheck.setElem(6, 1, 1);
	matrixCheck.setElem(6, 2, 1);
	matrixCheck.setElem(6, 3, 1);
	matrixCheck.setElem(7, 0, 1);
	matrixCheck.setElem(7, 1, 1);
	matrixCheck.setElem(7, 2, 1);
	matrixCheck.setElem(7, 3, 1);
	matrixCheck.setElem(8, 0, 1);
	matrixCheck.setElem(8, 1, 1);
	matrixCheck.setElem(8, 2, 1);
	matrixCheck.setElem(8, 3, 1);

	matrixCheck.setElem(4, 5, 1);
	matrixCheck.setElem(5, 4, 1);
	matrixCheck.setElem(6, 7, 1);
	matrixCheck.setElem(7, 6, 1);
	matrixCheck.setElem(7, 8, 1);
	matrixCheck.setElem(8, 7, 1);

	graph.set_numNode(9);
	graph.set_numLine(9);
	graph.set_numColumn(9);
	graph.set_Matrix(matrixCheck);
	matrixCheck.PrintMatrix();*/
	

	/*
	Matrix matrixCheck(numLine + 4, numColumn + 4);
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			if (i == j) {
				matrixCheck.setElem(i, j, 0);
				//arr[numColumn * i + j] = 0;
			}
			else {
				matrixCheck.setElem(i, j, 1);
				//arr[numColumn * i + j] = 1;
			}
		}
	}
	matrixCheck.setElem(0, 3, 1);
	matrixCheck.setElem(0, 4, 1);
	matrixCheck.setElem(0, 5, 1);
	matrixCheck.setElem(0, 6, 1);

	matrixCheck.setElem(1, 3, 1);
	matrixCheck.setElem(1, 4, 1);
	matrixCheck.setElem(1, 5, 1);
	matrixCheck.setElem(1, 6, 1);

	matrixCheck.setElem(2, 3, 1);
	matrixCheck.setElem(2, 4, 1);
	matrixCheck.setElem(2, 5, 1);
	matrixCheck.setElem(2, 6, 1);

	matrixCheck.setElem(3, 0, 1);
	matrixCheck.setElem(3, 1, 1);
	matrixCheck.setElem(3, 2, 1);
	matrixCheck.setElem(3, 4, 1);

	matrixCheck.setElem(4, 0, 1);
	matrixCheck.setElem(4, 1, 1);
	matrixCheck.setElem(4, 2, 1);
	matrixCheck.setElem(4, 3, 1);

	matrixCheck.setElem(5, 0, 1);
	matrixCheck.setElem(5, 1, 1);
	matrixCheck.setElem(5, 2, 1);
	matrixCheck.setElem(5, 6, 1);

	matrixCheck.setElem(6, 0, 1);
	matrixCheck.setElem(6, 1, 1);
	matrixCheck.setElem(6, 2, 1);
	matrixCheck.setElem(6, 5, 1);

	graph.set_numNode(7);
	graph.set_numLine(7);
	graph.set_numColumn(7);
	graph.set_Matrix(matrixCheck);
	matrixCheck.PrintMatrix();*/
	

	std::cout << std::endl;
	res = graph.algorithmEven();
	if (res) {
		arr_graph.push_back(graph);
	}
	fun(graph, &num, k, &arr_graph);
	int numVector = 0;

	//check_k_connected(arr_graph);	
	
	std::cout << std::endl;
	bool r = false;
	std::string str = "";
	std::vector<Graph> min_graph;
	for (auto v : arr_graph) {
		r = false;
		r = v.checkMinGraph();
		std::cout << std::endl;
		if (r == true) {
			v.set_minimality(r);
			min_graph.push_back(v);
		}
	}
	
	std::vector<Graph> finish_graph;
	for (auto m_graph : min_graph) {
		r = false;
		//m_graph.get_Matrix().PrintMatrix();
		r = m_graph.checkContractionMinmality();
		std::cout << std::endl;
		if (r == true) {
			m_graph.set_minContraction(r);
			finish_graph.push_back(m_graph);
		}
	}
	for (auto f_graph : finish_graph) {
		std::cout << "Matrix number: " << numVector << std::endl;
		Matrix m = f_graph.get_Matrix();
		m.PrintMatrix();
		int numEdge = m.currentNumEdge();
		std::cout << "numEdge " << numEdge << std::endl;
		std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
		std::cout << "Result of check min graph = true" << std::endl;
		std::cout << "Result of check min contraction graph = true" << std::endl;
		std::cout << std::endl;
		numVector++;
	}
	std::cout << "Num " << num + 1 << std::endl;
	std::cout << std::endl;
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	printToFile(finish_graph, "out.txt");
	arr_graph.clear();
	min_graph.clear();
	finish_graph.clear();
	check();
	return 0;
}
