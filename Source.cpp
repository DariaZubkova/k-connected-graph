#include "matrix.h"
#include <map>
#include <utility>

void extraMatrix(Matrix &matrix, int numI, int numJ) {
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	int extraVertex = 0;
	std::map<int, int> n;
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
	//std::vector<int> newArr(newNumLine * newNumColumn);
	for (int i = 0; i < newNumLine; i++) {
		for (int j = 0; j < newNumColumn; j++) {
			if (i != numI && i != numJ) {
				if (i < numLine) {
					//map
					auto prevJ = n.find(i);
					if (prevJ != n.end()) {
						if (j == prevJ->second) {//map
							newMatrix.setElem(i, j, 1);
						}
						else {
							newMatrix.setElem(i, j, 0);
						}
					}
				}
				else {
					if (j < numColumn) {
						//map
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
					int value = matrix.getElem(j, j);
					newMatrix.setElem(i, j, value);
					//newArr[newNumColumn * i + j] = arr[numColumn * i + j];
				}
				else {
					newMatrix.setElem(i, j, 0);
				}
			}
		}
	}
	newMatrix.PrintMatrix();
	std::cout << std::endl;
	n.clear();
}

void algorithmEven(Matrix &matrix) {
	extraMatrix(matrix, 0, 2);
}

void fun(Matrix &matrix, int* num) {
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				algorithmEven(matrix);
				fun(matrix, num);
				//check
				//algorithmEven(arr, numLine, numColumn);
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				//std::cout << "Return " << i << " " << j << std::endl;
				//PrintMatrix(arr, numLine, numColumn);
			}
		}
	}
	int st = 1;
	for (int i = 0; i < numLine; i++) {
		for (int j = 0; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1) {
				st = 0;
			}
		}
	}
	(*num) += st;
}


int main(void) {
	int numNode = 3;
	int numEdge = 0;
	int numLine = numNode;
	int numColumn = numNode;
	//std::vector<int> arr(numLine * numColumn);
	Matrix matrix(numLine, numColumn);
	std::vector<node>edgeBack;
	std::vector<edge>edgeExit;
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
	matrix.PrintMatrix();
	std::cout << "numEdge " << numEdge << std::endl;
	algorithmEven(matrix);
	fun(matrix, &num);
	std::cout << "Num " << num << std::endl;
	return 0;
}
