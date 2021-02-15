#include "matrix.h"
#include"graph.h"
#include <map>
#include <queue>
#include <utility>

//const int INF = 1000000000;


bool compareMatrix(std::vector<Matrix>& arr_matr, Matrix& matrix) {
	int nL = matrix.getNumLine(), nC = matrix.getNumColumn();
	bool res = true;
	//matrix.PrintMatrix();
	//std::cout << "COMPARE:" << std::endl;
	for (auto v : arr_matr) {
		res = true;
		//v.PrintMatrix();
		//std::cout << std::endl;
		for (int i = 0; i < nL; i++) {
			for (int j = 0; j < nC; j++) {
				if (v.getElem(i, j) != matrix.getElem(i, j)) {
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

void fun(Graph& graph, int* num, int k, std::vector<Matrix>* arr_matr) { //������� ������
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
				if (res && !compareMatrix((*arr_matr), matrix)) {
					std::cout << *num + 2 << std::endl;
					(*arr_matr).push_back(matrix);
				}
				(*num) += 1;
				fun(graph, num, k, arr_matr);
				//check
				//algorithmEven(arr, numLine, numColumn);
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "Return " << i << " " << j << std::endl;
				//PrintMatrix(arr, numLine, numColumn);
			}
		}
	}
}

int main(void) {
	int numNode = 6;
	int numEdge = 0;
	int k = 4;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);
	//std::vector<int> arr(numLine * numColumn);
	/*Matrix matrix(numLine, numColumn);
	Matrix f(numLine, numColumn);
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			f.setElem(i, j, 0);
		}
	}
	std::vector<int> ptr(numNode);
	std::vector<int> d(numNode);
	int maxFlow = 0;
	int num = 0;
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			if (i == j) {
				matrix.setElem(i, j, 0);
				//arr[numColumn * i + j] = 0;
			}
			else {
				matrix.setElem(i, j, 1);
				//arr[numColumn * i + j] = 1;
			}
		}
	}
	numEdge = matrix.currentNumEdge();
	std::cout << "Matrix: " << std::endl;
	matrix.PrintMatrix();
	std::cout << "numEdge " << numEdge << std::endl;
	std::cout << "Start algorithm " << std::endl;
	std::cout << std::endl;
	//maxFlow = dinic(matrix, f, ptr, d);
	
	Matrix matrixCheck(numLine, numColumn);
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			matrixCheck.setElem(i, j, 0);
		}
	}
	matrixCheck.setElem(0, 1, 1);
	matrixCheck.setElem(0, 2, 1);
	matrixCheck.setElem(0, 3, 1);
	matrixCheck.setElem(1, 0, 1);
	matrixCheck.setElem(1, 3, 1);
	matrixCheck.setElem(2, 0, 1);
	matrixCheck.setElem(3, 0, 1);
	matrixCheck.setElem(3, 1, 1);

	bool res;
	std::vector<Matrix> arr_matr;
	//algorithmEven(matrix, -1, k);
	res = algorithmEven(matrix, -1, k);
	if (res) {
		arr_matr.push_back(matrix);
	}*/
	//std::cout << std::endl;

	std::cout << "Start algorithm " << std::endl;
	std::cout << std::endl;
	bool res;
	int num = 0;
	std::vector<Matrix> arr_matr;
	res = graph.algorithmEven();
	if (res) {
		arr_matr.push_back(graph.get_Matrix());
	}
	fun(graph, &num, k, &arr_matr);
	int numVector = 0;
	//std::vector<Matrix> arr_matr;

	/*graph.algorithmEven(graph);
	int num = 0;
	Matrix matrix = graph.get_Matrix();
	fun(matrix, &num, k, &arr_matr);*/
	std::cout << std::endl;
	for (auto v : arr_matr) {
		std::cout << "Matrix number: " << numVector << std::endl;
		v.PrintMatrix();
		int numEdge = v.currentNumEdge();
		std::cout << "numEdge " << numEdge << std::endl;
		std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
		std::cout << std::endl;
		numVector++;
	}
	std::cout << "Num " << num + 1 << std::endl;
	//std::cout << "maxFlow " << maxFlow << std::endl;
	return 0;
}
