#pragma once
#include <vector>
#include<iostream>

typedef struct node {
	int i;
	int j;
}node;

typedef struct edge {
	int id;
	int i;
	int j;
}edge;

class Matrix {
public:
	Matrix(int numLine, int numColumn);
	Matrix() {
		numLine = 0;
		numColumn = 0;
	} //инициализация
	void printMatrix();
	int currentNumEdge();
	int getElem(int i, int j);
	void setElem(int i, int j, int value);
	void addVal(int i, int j, int value);
	void difVal(int i, int j, int value);
	int getNumLine();
	int getNumColumn();
	void setNumLine(int num);
	void setNumColumn(int num);
	void nullArray();
	~Matrix();
private:
	std::vector<int> arr;
	int numLine;
	int numColumn;
};