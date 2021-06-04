#include"graph.h"
#include"matrix.h"
#include <ctime>
#include <fstream>
#include <thread>
#include <mutex>

std::mutex mtx;
std::vector<Graph> array_graph;

bool compareMatrix(std::vector<Graph>& arr_graph, Matrix& matrix) {
	int nL = matrix.getNumLine(), nC = matrix.getNumColumn();
	bool res = true;
	if (arr_graph.empty()) {
		return false;
	}
	for (int i = 0; i < arr_graph.size(); i++) { // auto v : arr_graph
		Graph v = arr_graph[i];
		res = true;
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
	return res;
}

void fun_fun2(Graph& graph, int* num, int k, std::vector<Graph>* arr_graph) { //ïåðåáîð ãðàôîâ
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
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
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				dVertex[j]--;
				graph.set_Matrix(matrix);
				res = graph.algorithmEven();
				if (res) {
					(*arr_graph).push_back(graph);
				}
				(*num) += 1;
				fun_fun2(graph, num, k, arr_graph);
				dVertex[i]++;
				dVertex[j]++;
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}

void fun(Graph& graph, int* num, int k, std::vector<Graph>* arr_graph) { //ïåðåáîð ãðàôîâ
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<int> dVertex;
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex.push_back(s);
	}
	for (int i = 0; i < 1; i++) { //numLine - 1
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				dVertex[i]--;
				dVertex[j]--;
				graph.set_Matrix(matrix);
				res = graph.algorithmEven();
				if (res && !compareMatrix((*arr_graph), matrix)) {
					(*arr_graph).push_back(graph);
				}
				(*num) += 1;
				fun_fun2(graph, num, k, arr_graph);
				dVertex[i]++;
				dVertex[j]++;
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
}

void f(Graph graph, int k) {
	std::cout << "ID THREAD = " << std::this_thread::get_id() << std::endl;
	//if (graph.matrix.arr.empty()) {
	//	std::cout << "GOVNO =  " << std::this_thread::get_id() << std::endl;
	//	return;
	//}
	std::vector<std::thread> threads;
	//std::cout << "matrix get" << std::endl;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	//std::cout << "dVertex init" << std::endl;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine - 1; i++) {
		std::cout << "For i = " << i << std::endl;
		for (int j = i + 1; j < numColumn; j++) {
			std::cout << "For j = " << j << std::endl;
			//std::cout << "If" << std::endl;
			//matrix.PrintMatrix();
			//for (int re = 0; re < numLine; re++) {
			//	std::cout << "dVertex[" << re << "] = " << dVertex[re] << std::endl;
			//}
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) { // problem!
				//std::cout << "sets" << std::endl;
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				//std::cout << "i = " << i << std::endl;
				//std::cout << "j = " << j << std::endl;
				//std::cout << "set matrix" << std::endl;
				graph.set_Matrix(matrix);
				//std::cout << "set matrix end" << std::endl;
				//(*num)++;
				//std::cout << "lock" << std::endl;
				//std::lock_guard<std::mutex> lock(mtx);
				//std::cout << "Even" << std::endl;
				res = graph.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, matrix)) {
					//std::cout << *num + 1 << std::endl;
					//std::lock_guard<std::mutex> lock(mtx);
					//mutex_graph.lock();
					//std::cout << "push res" << std::endl;
					array_graph.push_back(graph);
					//mutex_graph.unlock();
				}

				//std::cout << "Next = " << std::this_thread::get_id() << std::endl;
				//std::cout << "recursive" << std::endl;
				f(graph, k);

				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
				//std::cout << "End = " << i << " " << j << std::endl;
			}
		}
	}
	dVertex.clear();
}

void thread_graph(Graph& graph, int k) {
	std::cout << "ID THREAD MAIN" << std::endl;
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<std::thread> threads;
	std::vector<int> dVertex(numLine);
	for (int i = 0; i < numLine - 1; i++) { //numLine - 1
		for (int j = i + 1; j < numColumn; j++) { //numColumn
			Graph gr(graph.get_numNode(), graph.get_k());
			Matrix m(numLine, numColumn);
			for (int u = 0; u < numLine; u++) {
				for (int l = 0; l < numColumn; l++) {
					m.setElem(u, l, matrix.getElem(u, l));
				}
			}
			//Matrix m = graph.get_Matrix();
			gr.set_Matrix(m);
			for (int h = 0; h < numLine; h++) {
				int s = 0;
				for (int l = 0; l < numColumn; l++) {
					s += matrix.getElem(h, l);
				}
				dVertex[h] = s;
			}
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				m.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
					m.setElem(j, i, 0);
				}
				graph.set_Matrix(matrix);
				gr.set_Matrix(m);
				//std::cout << "Delete " << i << " " << j << std::endl;
				//(*num)++;
				//PrintMatrix(arr, numLine, numColumn);
				//check
				//res = graph.algorithmEven();
				res = gr.algorithmEven();
				//res = algo(graph);
				if (res && !compareMatrix(array_graph, m)) { //matrix
					//std::cout << *num + 1 << std::endl;
					std::lock_guard<std::mutex> lock(mtx);
					array_graph.push_back(gr); //graph

				}

				std::thread thr(f, gr, k); //graph
				//threads.push_back(std::move(thr));
				threads.emplace_back(std::move(thr));
				//threads.emplace_back(fun_thread, gr, k);
				//std::thread t(&fun_thread, std::ref(gr), k);
				//t.detach();
				//threads_graph.push_back(std::move(t));
				//threads_graph.push_back(std::thread(&fun_thread, std::ref(gr), k, arr_graph));
				//fun_thread(gr, k, arr_graph);

				matrix.setElem(i, j, 1);
				m.setElem(i, j, 0);
				matrix.setElem(j, i, 1);
				m.setElem(j, i, 0);
				graph.set_Matrix(matrix);
				gr.set_Matrix(m);
			}
			//dVertex.clear();
		}
	}
	//dVertex.clear();
	for (auto& thr : threads) {
		thr.join();
	}
	//for (int i = 0; i < threads.size(); i++) {
	//	threads_graph[i].detach();
	//}
}

void enumeration_graph(std::vector<Graph> all_arr_graph, int k) {
	while (!all_arr_graph.empty()) {
		Graph graph = all_arr_graph[0];
		all_arr_graph.erase(all_arr_graph.begin());
		Matrix matrix = graph.get_Matrix();
		int numLine = matrix.getNumLine();
		int numColumn = matrix.getNumColumn();
		if (numLine > 0 || numColumn > 0) {
			bool res = false;
			std::vector<int> dVertex(numLine, 0);
			for (int h = 0; h < numLine; h++) {
				int s = 0;
				for (int l = 0; l < numColumn; l++) {
					s += matrix.getElem(h, l);
				}
				dVertex[h] = s;
			}
			for (int i = 0; i < numLine - 1; i++) {
				for (int j = i + 1; j < numColumn; j++) {
					if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
						matrix.setElem(i, j, 0);
						if (matrix.getElem(j, i) == 1) {
							matrix.setElem(j, i, 0);
						}
						graph.set_Matrix(matrix);
						//res = graph.algorithmGalil();
						res = graph.algorithmEven();
						if (res) {
							std::lock_guard<std::mutex> lock(mtx);
							array_graph.push_back(graph);

						}
						all_arr_graph.push_back(graph);
						matrix.setElem(i, j, 1);
						matrix.setElem(j, i, 1);
						graph.set_Matrix(matrix);
					}
				}
			}
			dVertex.clear();
		}
	}
}

void start_enumeration_graph(std::vector<Graph> all_arr_graph, int k) {
	Graph graph = all_arr_graph[0];
	all_arr_graph.erase(all_arr_graph.begin());
	Matrix matrix = graph.get_Matrix();
	int numLine = matrix.getNumLine();
	int numColumn = matrix.getNumColumn();
	bool res = false;
	std::vector<std::thread> threads;
	std::vector<int> dVertex(numLine, 0);
	for (int h = 0; h < numLine; h++) {
		int s = 0;
		for (int l = 0; l < numColumn; l++) {
			s += matrix.getElem(h, l);
		}
		dVertex[h] = s;
	}
	for (int i = 0; i < 1; i++) { //numLine - 1
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1 && dVertex[i] > k && dVertex[j] > k) {
				matrix.setElem(i, j, 0);
				if (matrix.getElem(j, i) == 1) {
					matrix.setElem(j, i, 0);
				}
				graph.set_Matrix(matrix);
				//res = graph.algorithmGalil();
				res = graph.algorithmEven();
				if (res) {
					std::lock_guard<std::mutex> lock(mtx);
					array_graph.push_back(graph);

				}
				std::vector<Graph> graph_func;
				graph_func.push_back(graph);
				std::thread thr(enumeration_graph, graph_func, k);
				threads.emplace_back(std::move(thr));
				//perebor(graph_func, k);
				graph_func.erase(graph_func.begin());
				graph_func.clear();
				matrix.setElem(i, j, 1);
				matrix.setElem(j, i, 1);
				graph.set_Matrix(matrix);
			}
		}
	}
	dVertex.clear();
	for (auto& thr : threads) {
		thr.join();
	}
}

void check_k_connected(std::vector<Graph> arr) {
	int numVector = 0;
	if (!arr.empty()) {
		int k = arr[0].get_k();
		for (auto v : arr) {
			std::cout << "Matrix number: " << numVector << std::endl;
			Matrix m = v.get_Matrix();
			m.printMatrix();
			int numEdge = m.currentNumEdge();
			std::cout << "numEdge " << numEdge << std::endl;
			std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
			std::cout << std::endl;
			numVector++;
		}
	}
	std::cout << "END!" << std::endl;
}

void fromFileInGraph(std::string nameFile) {
	std::fstream f;
	Graph graph(20, 5);
	Matrix matrix(20, 20);
	f.open(nameFile, std::ios::out);
	if (f)
	{
		std::string str;
		int numNode = 0;
		int flag = 0;
		int i = 0, j = 0, value = 1;
		while (std::getline(f, str)) {
			if (str[0] == 'T') {
				std::string s = "";
				int l = 2;
				while (str[l] != '\n') {
					s += str[l];
				}
				numNode = std::stoi(s) + 1;
			}
			if (str[0] == 'E') {
				if (flag == 0) {
					; //make matrix
				}
				flag = 1;
				int l = 2;
			}
			
		}
		f.close();
	}
}

void printToFile(std::vector<Graph> graphs, std::string nameFile) {
	std::fstream f;
	f.open(nameFile, std::ios::out);
	if (f)
	{
		std::string str;
		int number = 1;
		for (auto f_graph : graphs)
		{
			int numLine = f_graph.get_numLine();
			int numColumn = f_graph.get_numColumn();
			Matrix m = f_graph.get_Matrix();
			str = "N ";
			str += std::to_string(number);
			f << str << std::endl;
			str = "";
			for (int i = 0; i < numLine; i++) {
				str = "T ";
				str += std::to_string(i);
				f << str << std::endl;
				str = "";
			}
			for (int i = 0; i < numLine; i++) {
				for (int j = i + 1; j < numColumn; j++) {
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
			number++;
			f << std::endl;
		}
		f.close();
	}
}

void check0() {
	int numNode = 15;
	int numEdge = 0;
	int k = 4;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);

	Matrix matr(numLine, numColumn);

	//0
	matr.setElem(0, 1, 1);
	matr.setElem(1, 0, 1);
	matr.setElem(0, 2, 1);
	matr.setElem(2, 0, 1);
	matr.setElem(0, 3, 1);
	matr.setElem(3, 0, 1);
	matr.setElem(0, 4, 1);
	matr.setElem(4, 0, 1);
	//1
	matr.setElem(1, 2, 1);
	matr.setElem(2, 1, 1);
	matr.setElem(1, 5, 1);
	matr.setElem(5, 1, 1);
	matr.setElem(1, 6, 1);
	matr.setElem(6, 1, 1);
	//3
	matr.setElem(2, 5, 1);
	matr.setElem(5, 2, 1);
	matr.setElem(2, 6, 1);
	matr.setElem(6, 2, 1);
	//4
	matr.setElem(3, 7, 1);
	matr.setElem(7, 3, 1);
	matr.setElem(3, 8, 1);
	matr.setElem(8, 3, 1);
	matr.setElem(3, 11, 1);
	matr.setElem(11, 3, 1);
	matr.setElem(3, 13, 1);
	matr.setElem(13, 3, 1);
	//5
	matr.setElem(4, 7, 1);
	matr.setElem(7, 4, 1);
	matr.setElem(4, 8, 1);
	matr.setElem(8, 4, 1);
	matr.setElem(4, 11, 1);
	matr.setElem(11, 4, 1);
	matr.setElem(4, 13, 1);
	matr.setElem(13, 4, 1);
	//7
	matr.setElem(5, 9, 1);
	matr.setElem(9, 5, 1);
	matr.setElem(5, 10, 1);
	matr.setElem(10, 5, 1);
	matr.setElem(5, 12, 1);
	matr.setElem(12, 5, 1);
	matr.setElem(5, 14, 1);
	matr.setElem(14, 5, 1);
	//8
	matr.setElem(6, 9, 1);
	matr.setElem(9, 6, 1);
	matr.setElem(6, 10, 1);
	matr.setElem(10, 6, 1);
	matr.setElem(6, 12, 1);
	matr.setElem(12, 6, 1);
	matr.setElem(6, 14, 1);
	matr.setElem(14, 6, 1);
	//9
	matr.setElem(7, 8, 1);
	matr.setElem(8, 7, 1);
	matr.setElem(7, 9, 1);
	matr.setElem(9, 7, 1);
	matr.setElem(7, 10, 1);
	matr.setElem(10, 7, 1);
	//10
	matr.setElem(8, 10, 1);
	matr.setElem(10, 8, 1);
	//11
	matr.setElem(9, 10, 1);
	matr.setElem(10, 9, 1);
	//12
	//13
	matr.setElem(11, 13, 1);
	matr.setElem(13, 11, 1);
	matr.setElem(11, 14, 1);
	matr.setElem(14, 11, 1);
	//14
	matr.setElem(12, 13, 1);
	matr.setElem(13, 12, 1);
	matr.setElem(12, 14, 1);
	matr.setElem(14, 12, 1);
	//15
	matr.setElem(13, 14, 1);
	matr.setElem(14, 13, 1);
	//16

	graph.set_Matrix(matr);
	std::vector<Graph> gr;
	gr.push_back(graph);
	printToFile(gr, "out_gr.txt");

	unsigned int start_time = clock();
	bool res = false;
	res = graph.algorithmEven();
	//res = graph.algorithmGalil();
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	if (res) {
		std::cout << "Result k-coonected = true" << std::endl;
	}
	else {
		std::cout << "Result k-coonected = false" << std::endl;
	}
	std::cout << "Time work = " << search_time << std::endl;
	//std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	res = graph.checkMinGraph();
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - end_time;
	if (res) {
		std::cout << "Result min = true" << std::endl;
	}
	else {
		std::cout << "Result min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time2 << std::endl;
	//std::cout << "Time work = " << search_time2 / CLOCKS_PER_SEC << std::endl;
	res = graph.checkContractionMinmality();
	//res = graph.perebor_vertex_check();
	unsigned int end_time3 = clock();
	unsigned int search_time3 = end_time3 - end_time2;
	if (res) {
		std::cout << "Result contr min = true" << std::endl;
	}
	else {
		std::cout << "Result contr min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time3 << std::endl;
	//std::cout << "Time work = " << search_time3 / CLOCKS_PER_SEC << std::endl;
}

void check() {
	int numNode = 17;
	int numEdge = 0;
	int k = 5;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);

	Matrix matr(numLine, numColumn);
	/*for (int i = 0; i < numLine; i++) {
		for (int j = i+6; j < numColumn; j+=6) {
			matr.setElem(i, j, 1);
			matr.setElem(j, i, 1);
		}
		if (i < (numNode - 1)) {
			matr.setElem(i, i + 1, 1);
			matr.setElem(i + 1, i, 1);
		}
	}
	matr.setElem(0, numNode - 1, 1);
	matr.setElem(numNode - 1, 0, 1);*/

	//0
	for (int i = 1; i < 6; i++) {
		matr.setElem(0, i, 1);
		matr.setElem(i, 0, 1);
	}
	//1
	matr.setElem(1, 2, 1);
	matr.setElem(2, 1, 1);
	matr.setElem(1, 3, 1);
	matr.setElem(3, 1, 1);
	matr.setElem(1, 7, 1);
	matr.setElem(7, 1, 1);
	matr.setElem(1, 8, 1);
	matr.setElem(8, 1, 1);
	//2
	matr.setElem(2, 4, 1);
	matr.setElem(4, 2, 1);
	matr.setElem(2, 5, 1);
	matr.setElem(5, 2, 1);
	matr.setElem(2, 6, 1);
	matr.setElem(6, 2, 1);
	//3
	matr.setElem(3, 6, 1);
	matr.setElem(6, 3, 1);
	matr.setElem(3, 7, 1);
	matr.setElem(7, 3, 1);
	matr.setElem(3, 8, 1);
	matr.setElem(8, 3, 1);
	//4
	matr.setElem(4, 9, 1);
	matr.setElem(9, 4, 1);
	matr.setElem(4, 10, 1);
	matr.setElem(10, 4, 1);
	matr.setElem(4, 13, 1);
	matr.setElem(13, 4, 1);
	matr.setElem(4, 15, 1);
	matr.setElem(15, 4, 1);
	//5
	matr.setElem(5, 9, 1);
	matr.setElem(9, 5, 1);
	matr.setElem(5, 10, 1);
	matr.setElem(10, 5, 1);
	matr.setElem(5, 13, 1);
	matr.setElem(13, 5, 1);
	matr.setElem(5, 15, 1);
	matr.setElem(15, 5, 1);
	//6
	matr.setElem(6, 10, 1);
	matr.setElem(10, 6, 1);
	matr.setElem(6, 11, 1);
	matr.setElem(11, 6, 1);
	matr.setElem(6, 13, 1);
	matr.setElem(13, 6, 1);
	matr.setElem(6, 14, 1);
	matr.setElem(14, 6, 1);
	//7
	matr.setElem(7, 11, 1);
	matr.setElem(11, 7, 1);
	matr.setElem(7, 12, 1);
	matr.setElem(12, 7, 1);
	matr.setElem(7, 14, 1);
	matr.setElem(14, 7, 1);
	matr.setElem(7, 16, 1);
	matr.setElem(16, 7, 1);
	//8
	matr.setElem(8, 11, 1);
	matr.setElem(11, 8, 1);
	matr.setElem(8, 12, 1);
	matr.setElem(12, 8, 1);
	matr.setElem(8, 14, 1);
	matr.setElem(14, 8, 1);
	matr.setElem(8, 16, 1);
	matr.setElem(16, 8, 1);
	//9
	matr.setElem(9, 10, 1);
	matr.setElem(10, 9, 1);
	matr.setElem(9, 11, 1);
	matr.setElem(11, 9, 1);
	matr.setElem(9, 12, 1);
	matr.setElem(12, 9, 1);
	//10
	matr.setElem(10, 12, 1);
	matr.setElem(12, 10, 1);
	//11
	matr.setElem(11, 12, 1);
	matr.setElem(12, 11, 1);
	//12
	//13
	matr.setElem(13, 15, 1);
	matr.setElem(15, 13, 1);
	matr.setElem(13, 16, 1);
	matr.setElem(16, 13, 1);
	//14
	matr.setElem(14, 15, 1);
	matr.setElem(15, 14, 1);
	matr.setElem(14, 16, 1);
	matr.setElem(16, 14, 1);
	//15
	matr.setElem(15, 13, 1);
	matr.setElem(13, 15, 1);
	matr.setElem(15, 14, 1);
	matr.setElem(14, 15, 1);
	matr.setElem(15, 16, 1);
	matr.setElem(16, 15, 1);
	//16

	graph.set_Matrix(matr);
	std::vector<Graph> gr;
	gr.push_back(graph);
	printToFile(gr, "out_gr.txt");

	unsigned int start_time = clock();
	bool res = false;
	//res = graph.algorithmEven();
	res = graph.algorithmGalil();
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	if (res) {
		std::cout << "Result k-coonected = true" << std::endl;
	}
	else {
		std::cout << "Result k-coonected = false" << std::endl;
	}
	std::cout << "Time work = " << search_time << std::endl;
	//std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	res = graph.checkMinGraph();
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - end_time;
	if (res) {
		std::cout << "Result min = true" << std::endl;
	}
	else {
		std::cout << "Result min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time2 << std::endl;
	//std::cout << "Time work = " << search_time2 / CLOCKS_PER_SEC << std::endl;
	res = graph.checkContractionMinmality();
	//res = graph.perebor_vertex_check();
	unsigned int end_time3 = clock();
	unsigned int search_time3 = end_time3 - end_time2;
	if (res) {
		std::cout << "Result contr min = true" << std::endl;
	}
	else {
		std::cout << "Result contr min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time3 << std::endl;
	//std::cout << "Time work = " << search_time3 / CLOCKS_PER_SEC << std::endl;
}

void check4() {
	int numNode = 6;
	int numEdge = 0;
	int k = 4;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);

	Matrix matr(numLine, numColumn);
	//0
	matr.setElem(0, 1, 1);
	matr.setElem(1, 0, 1);
	matr.setElem(0, 2, 1);
	matr.setElem(2, 0, 1);
	matr.setElem(0, 4, 1);
	matr.setElem(4, 0, 1);
	matr.setElem(0, 5, 1);
	matr.setElem(5, 0, 1);
	//1
	matr.setElem(1, 2, 1);
	matr.setElem(2, 1, 1);
	matr.setElem(1, 3, 1);
	matr.setElem(3, 1, 1);
	matr.setElem(1, 5, 1);
	matr.setElem(5, 1, 1);
	//2
	matr.setElem(2, 3, 1);
	matr.setElem(3, 2, 1);
	matr.setElem(2, 4, 1);
	matr.setElem(4, 2, 1);
	//3
	matr.setElem(3, 4, 1);
	matr.setElem(4, 3, 1);
	matr.setElem(3, 5, 1);
	matr.setElem(5, 3, 1);
	//4
	matr.setElem(4, 5, 1);
	matr.setElem(5, 4, 1);

	graph.set_Matrix(matr);
	std::vector<Graph> gr;
	gr.push_back(graph);
	printToFile(gr, "out_gr.txt");


	unsigned int start_time = clock();
	bool res = false;

	res = graph.algorithmEven();
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	if (res) {
		std::cout << "Result k-coonected = true" << std::endl;
	}
	else {
		std::cout << "Result k-coonected = false" << std::endl;
	}
	std::cout << "Time work = " << search_time << std::endl;
	//std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	res = graph.checkMinGraph();
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - end_time;
	if (res) {
		std::cout << "Result min = true" << std::endl;
	}
	else {
		std::cout << "Result min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time2 << std::endl;
	//std::cout << "Time work = " << search_time2 / CLOCKS_PER_SEC << std::endl;
	res = graph.checkContractionMinmality();
	unsigned int end_time3 = clock();
	unsigned int search_time3 = end_time3 - end_time2;
	if (res) {
		std::cout << "Result contr min = true" << std::endl;
	}
	else {
		std::cout << "Result contr min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time3 << std::endl;
	//std::cout << "Time work = " << search_time3 / CLOCKS_PER_SEC << std::endl;
}

void check77() {
	int numNode = 15;
	int numEdge = 0;
	int k = 6;
	int numLine = numNode;
	int numColumn = numNode;
	Graph graph(numNode, k);

	Matrix matr(numLine, numColumn);

	//0
	matr.setElem(0, 1, 1);
	matr.setElem(1, 0, 1);
	matr.setElem(0, 2, 1);
	matr.setElem(2, 0, 1);
	matr.setElem(0, 3, 1);
	matr.setElem(3, 0, 1);
	matr.setElem(0, 4, 1);
	matr.setElem(4, 0, 1);
	matr.setElem(0, 7, 1);
	matr.setElem(7, 0, 1);
	matr.setElem(0, 8, 1);
	matr.setElem(8, 0, 1);
	//1
	matr.setElem(1, 2, 1);
	matr.setElem(2, 1, 1);
	matr.setElem(1, 3, 1);
	matr.setElem(3, 1, 1);
	matr.setElem(1, 4, 1);
	matr.setElem(4, 1, 1);
	matr.setElem(1, 5, 1);
	matr.setElem(5, 1, 1);
	matr.setElem(1, 6, 1);
	matr.setElem(6, 1, 1);
	//2
	matr.setElem(2, 5, 1);
	matr.setElem(5, 2, 1);
	matr.setElem(2, 6, 1);
	matr.setElem(6, 2, 1);
	matr.setElem(2, 7, 1);
	matr.setElem(7, 2, 1);
	matr.setElem(2, 8, 1);
	matr.setElem(8, 2, 1);
	//3
	matr.setElem(3, 9, 1);
	matr.setElem(9, 3, 1);
	matr.setElem(3, 11, 1);
	matr.setElem(11, 3, 1);
	matr.setElem(3, 12, 1);
	matr.setElem(12, 3, 1);
	matr.setElem(3, 14, 1);
	matr.setElem(14, 3, 1);
	//4
	matr.setElem(4, 9, 1);
	matr.setElem(9, 4, 1);
	matr.setElem(4, 11, 1);
	matr.setElem(11, 4, 1);
	matr.setElem(4, 12, 1);
	matr.setElem(12, 4, 1);
	matr.setElem(4, 14, 1);
	matr.setElem(14, 4, 1);
	//5
	matr.setElem(5, 9, 1);
	matr.setElem(9, 5, 1);
	matr.setElem(5, 10, 1);
	matr.setElem(10, 5, 1);
	matr.setElem(5, 12, 1);
	matr.setElem(12, 5, 1);
	matr.setElem(5, 13, 1);
	matr.setElem(13, 5, 1);
	//6
	matr.setElem(6, 9, 1);
	matr.setElem(9, 6, 1);
	matr.setElem(6, 10, 1);
	matr.setElem(10, 6, 1);
	matr.setElem(6, 12, 1);
	matr.setElem(12, 6, 1);
	matr.setElem(6, 13, 1);
	matr.setElem(13, 6, 1);
	//7
	matr.setElem(7, 10, 1);
	matr.setElem(10, 7, 1);
	matr.setElem(7, 11, 1);
	matr.setElem(11, 7, 1);
	matr.setElem(7, 13, 1);
	matr.setElem(13, 7, 1);
	matr.setElem(7, 14, 1);
	matr.setElem(14, 7, 1);
	//8
	matr.setElem(8, 10, 1);
	matr.setElem(10, 8, 1);
	matr.setElem(8, 11, 1);
	matr.setElem(11, 8, 1);
	matr.setElem(8, 13, 1);
	matr.setElem(13, 8, 1);
	matr.setElem(8, 14, 1);
	matr.setElem(14, 8, 1);
	//9
	matr.setElem(9, 10, 1);
	matr.setElem(10, 9, 1);
	matr.setElem(9, 11, 1);
	matr.setElem(11, 9, 1);
	//10
	matr.setElem(10, 11, 1);
	matr.setElem(11, 10, 1);
	//11
	//12
	matr.setElem(12, 13, 1);
	matr.setElem(13, 12, 1);
	matr.setElem(12, 14, 1);
	matr.setElem(14, 12, 1);
	//13
	matr.setElem(13, 14, 1);
	matr.setElem(14, 13, 1);
	//14

	graph.set_Matrix(matr);
	std::vector<Graph> gr;
	gr.push_back(graph);
	printToFile(gr, "out_gr.txt");

	unsigned int start_time = clock();
	bool res = false;
	res = graph.algorithmEven();
	//res = graph.algorithmGalil();
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	if (res) {
		std::cout << "Result k-coonected = true" << std::endl;
	}
	else {
		std::cout << "Result k-coonected = false" << std::endl;
	}
	std::cout << "Time work = " << search_time << std::endl;
	//std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	res = graph.checkMinGraph();
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - end_time;
	if (res) {
		std::cout << "Result min = true" << std::endl;
	}
	else {
		std::cout << "Result min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time2 << std::endl;
	//std::cout << "Time work = " << search_time2 / CLOCKS_PER_SEC << std::endl;
	res = graph.checkContractionMinmality();
	//res = graph.perebor_vertex_check();
	unsigned int end_time3 = clock();
	unsigned int search_time3 = end_time3 - end_time2;
	if (res) {
		std::cout << "Result contr min = true" << std::endl;
	}
	else {
		std::cout << "Result contr min = false" << std::endl;
	}
	std::cout << "Time work = " << search_time3 << std::endl;
	//std::cout << "Time work = " << search_time3 / CLOCKS_PER_SEC << std::endl;
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
	std::vector<Graph> all_arr_graph;
	//check3();
	//check6();
	check77();
	//check();
	/*Matrix matr(numLine, numColumn);
	matr.setElem(0, 1, 1);
	matr.setElem(0, 2, 1);
	matr.setElem(0, 4, 1);
	matr.setElem(0, 5, 1);
	matr.setElem(1, 0, 1);
	matr.setElem(1, 2, 1);
	matr.setElem(1, 3, 1);
	matr.setElem(1, 5, 1);
	matr.setElem(2, 0, 1);
	matr.setElem(2, 1, 1);
	matr.setElem(2, 3, 1);
	matr.setElem(2, 4, 1);
	matr.setElem(3, 1, 1);
	matr.setElem(3, 2, 1);
	matr.setElem(3, 4, 1);
	matr.setElem(3, 5, 1);
	matr.setElem(4, 0, 1);
	matr.setElem(4, 2, 1);
	matr.setElem(4, 3, 1);
	matr.setElem(4, 5, 1);
	matr.setElem(5, 0, 1);
	matr.setElem(5, 1, 1);
	matr.setElem(5, 3, 1);
	matr.setElem(5, 4, 1);
	graph.set_Matrix(matr);*/

	std::cout << std::endl;
	
	res = graph.algorithmEven();
	//res = graph.algorithmGalil();
	if (res) {
		//arr_graph.push_back(graph);
		array_graph.push_back(graph);
	}
	all_arr_graph.push_back(graph);

	//проверка на к-связность в переборе графов

	//enumeration_graph(all_arr_graph, k);
	start_enumeration_graph(all_arr_graph, k);
	//fun_fun2(graph, &num, k, &array_graph);

	//thread_graph(graph, k);
	int numVector = 0;

	//check_k_connected(arr_graph);	
	
	//arr_graph = array_graph;
	std::cout << std::endl;

	//проверка на одинаковые графы, одинаковые не добавляются в arr_graph
	bool res_compare = false;
	while (!array_graph.empty()) {
		Graph v = array_graph[0];
		array_graph.erase(array_graph.begin());
		res_compare = false;
		Matrix m = v.get_Matrix();
		res_compare = compareMatrix(array_graph, m);
		if (res_compare == false) {
			arr_graph.push_back(v);
		}
	}
	
	//проверка на минимальность
	bool r = false;
	std::string str = "";
	std::vector<Graph> min_graph;
	for (auto v : arr_graph) {
		r = false;
		r = v.checkMinGraph();
		if (r == true) {
			v.set_minimality(r);
			min_graph.push_back(v);
		}
	}

	//проверка на минимальность по стягиванию
	std::vector<Graph> finish_graph;
	for (auto m_graph : min_graph) {
		r = false;
		//r = m_graph.perebor_vertex_check();
		r = m_graph.checkContractionMinmality();
		if (r == true) {
			m_graph.set_minContraction(r);
			finish_graph.push_back(m_graph);
		}
	}


	//напечатать полученные графы
	for (auto f_graph : finish_graph) {
		int index_graph = 0;
		std::cout << "Matrix number: " << numVector << std::endl;
		Matrix m = f_graph.get_Matrix();
		std::vector<int> dVertex(numLine, 0);
		for (int h = 0; h < numLine; h++) {
			int s = 0;
			for (int l = 0; l < numColumn; l++) {
				s += m.getElem(h, l);
			}
			dVertex[h] = s;
		}
		int d_v = 0;
		for (int i = 0; i < numLine; i++) {
			if (dVertex[i] == k) {
				d_v++;
			}
		}
		m.printMatrix();
		int numEdge = m.currentNumEdge();
		std::cout << "numEdge " << numEdge << std::endl;
		std::cout << "CONGRATULATE!!!! It is k-connected graph; k = " << k << std::endl;
		std::cout << "Result of check min graph = true" << std::endl;
		std::cout << "Result of check min contraction graph = true" << std::endl;
		std::cout << "Number of vertex with d_v(k) = " << d_v << std::endl;
		std::cout << std::endl;
		index_graph++;
		numVector++;
	}
	//std::cout << "Num " << num + 1 << std::endl;
	std::cout << std::endl;
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	//std::cout << "Time work = " << search_time << std::endl;
	std::cout << "Time work = " << search_time / CLOCKS_PER_SEC << std::endl;
	printToFile(finish_graph, "out1.txt");
	arr_graph.clear();
	min_graph.clear();
	finish_graph.clear();

	//check();
	return 0;
}