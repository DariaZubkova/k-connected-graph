#include"graph.h"
#include <mutex>
#include <fstream>
std::mutex mtx_dinic;
#include <queue>
//std::lock_guard<std::mutex> lock(mtx_dinic);


const int inf = 2000;

int Graph::DFS(std::list <std::pair <int, std::pair <int, int> > >& block, std::vector <std::vector <int> >& lay, int& n, int& s, int& t, int cost)
{
	if (s == t) //если мы пришли в сток
		return cost; //вернем стоимость
	else
	{
		if (cost != 0)
		{
			int maxstream = 0; //это поток
			for (int i = 0; i < n; i++)
				if (lay[s][i] > 0) //если в остаточной сети есть ребро 
				{
					int flow = DFS(block, lay, n, i, t, std::min(cost, lay[s][i])); //строим поток
					if (flow) //если мы дошли до стока, то начинаем выстраивать блокирующий поток
					{
						block.push_back(std::make_pair(s, std::make_pair(i, flow)));
						lay[s][i] -= flow;
						maxstream += flow;
						cost -= flow;
					}
				}
			return (maxstream);
		}
		else return 0;
	}
}

void Graph::build(std::vector <std::vector <int> >& lay, int& n, int& s, int& t) //по остаточной сети строит слоистую сеть
{
	//инициализируем слоистую сеть
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			lay[i][j] = 0;
	std::vector <int> level(n, -1); //уровни для вершин
	level[s] = 0; //уровень истока 0
	std::queue <int> queue; //для обхода в ширину
	queue.push(s);
	int v;
	while (!queue.empty())
	{
		v = queue.front();
		queue.pop();
		for (int i = 0; i < n; i++)
			if (extra_matrix.getElem(v,i) > 0 && level[i] == -1) //если вершина есть, но уровень для нее пока не задан
			{
				level[i] = level[v] + 1; //увеличиваем для нее уровень
				queue.push(i);
			}
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (level[j] == level[i] + 1) //если вершины оказались на одном и том же уровне
				lay[i][j] = extra_matrix.getElem(i, j);
}

int Graph::Dinitz(int& s, int& t)
{
	int n = extra_matrix.getNumLine();
	std::vector <std::vector <int> > lay(n, std::vector <int>(n)); //это наша остаточная сеть,просто выделим память
	std::list <std::pair <int, std::pair <int, int> > > block; //список ребер блокирующего потока
	int stream = 0; //поток остаточной сети
	int flow = 0; //поток по ребрам блокирующего потока
	do
	{
		build(lay, n, s, t); //выстраиваем остаточную сеть
		flow = DFS(block, lay, n, s, t, inf); //flow - значение потока по ребрам блокирующего потока, block - блокирующий поток
		while (!block.empty()) //пока есть ребра в блокирующем потоке
		{
			std::pair <int, std::pair <int, int> > v = block.front();
			int a = block.front().first;
			int b = block.front().second.first;
			int flow = block.front().second.second;
			extra_matrix.difVal(a, b, flow);
			extra_matrix.addVal(b, a, flow);
			block.pop_front();
		}
		stream += flow; //блокирующий поток закончился, прибавляем его поток к потоку остаточной сети
	} while (flow); //пока блокирующие потоки,выполняем алгоритм
	return stream; //возвращаем поток остаточной сети
}




void printToFile2(Graph f_graph, std::string nameFile) {
	std::fstream f;
	f.open(nameFile, std::ios::out);
	if (f)
	{
		std::string str;
		int number = 1;
		Matrix extra_matrix = f_graph.get_ExtraMatrix();
		int numLine = extra_matrix.getNumLine();
		int numColumn = extra_matrix.getNumColumn();
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
			for (int j = 0; j < numColumn; j++) {
				if (extra_matrix.getElem(i, j) > 0) {
					str = "E ";
					str += std::to_string(i);
					str += " ";
					str += std::to_string(j);
					str += " ";
					str += std::to_string(1);
					f << str << std::endl;
					str = "";
				}
			}
		}
		number++;
		f << std::endl;
		f.close();
	}
}



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
	/*int j = newNumColumn - 1;
	for (int i = 0; i < newNumLine; i++) {
		if (i < J) {
			newMatrix.setElem(i, j, 1);
		}
	}
	int i = newNumLine - 1;
	for (int j = 0; j < newNumColumn; j++) {
		if (j < J) {
			newMatrix.setElem(i, j, 1);
		}
	}*/
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

int Graph::dfs(int u, int minFlow, Matrix& f, std::vector<int> ptr, std::vector<int> d, int t) {
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
		//std::cout << value_f << " " << value_extra << std::endl;
		if (d[v] == d[u] + 1) { //&& value_extra != 0
			int delta = dfs(v, min(minFlow, value_extra - value_f), f, ptr, d, t);
			if (delta > 0) {
				f.addVal(u, v, delta);
				f.difVal(v, u, delta);
				extra_matrix.difVal(u, v, delta);
				extra_matrix.addVal(v, u, delta);
				return delta;
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
	std::vector<int> d(numLine, -1);
	while (bfs(f, d, s, t)) {
		for (int i = 0; i < numLine; i++)
			ptr[i] = 0;
		//std::lock_guard<std::mutex> lock(mtx_dinic);
		flow = dfs(s, INF, f, ptr, d, t);
		while (flow != 0) {
			maxFlow += flow;
			if (maxFlow >= k) {
				return maxFlow;
			}
			flow = dfs(s, INF, f, ptr, d, t);
		}
	}
	return maxFlow;
}

bool Graph::algorithmEven() {
	int maxFlow = -1;
	if (numLine <= 0 || numColumn <= 0) {
		return false;
	}
	/*std::vector<int> dVertex(numLine, 0);
	for (int i = 0; i < numLine; i++) {
		int s = 0;
		for (int j = 0; j < numColumn; j++) {
			s += matrix.getElem(i, j);
		}
		dVertex[i] = s;
	}
	for (int i = 0; i < numLine; i++) {
		if (dVertex[i] < k) {
			return false;
		}
	}*/
	for (int j = 1; j < k; j++) {
		for (int i = 0; i < j; i++) {
			(*this).extraMatrix(i, j);
			//maxFlow = (*this).Dinitz(i, j);
			maxFlow = (*this).dinic(i, j);
			//maxFlow = gr.dinic(f, ptr, d, i, j);
			if (maxFlow < k) {
				return false;
			}
		}
	}
	int pos_s = matrix.getNumColumn();
	Graph newGraph = extraGraph((*this), k);
	for (int j = k; j < matrix.getNumColumn(); j++) {
		newGraph.extraMatrix(j, pos_s);
		//maxFlow = newGraph.Dinitz(j, pos_s);
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
			maxFlow = (*this).Dinitz(i, kk);
			//maxFlow = (*this).dinic(i, j);
			maxFlowVec.push_back(maxFlow);
		}
		int minFlow = 1000000;
		for (int i = 0; i < maxFlowVec.size(); i++) {
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
			//Graph newGraph = extraGraph((*this), j);
			newGraph.extraMatrix(j, pos_s);
			//printToFile2(newGraph, "named.txt");
			//newGraph.set_ExtraMatrix(matrix);
			//Matrix extraMatrix = newGraph.get_ExtraMatrix();
			maxFlow = newGraph.Dinitz(j, pos_s);
			//maxFlow = newGraph.dinic(j, pos_s);
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
					/*Matrix f(extra_matrix.getNumLine(), extra_matrix.getNumColumn());
					for (int i = 0; i < extra_matrix.getNumLine(); i++) {
						for (int j = 0; j < extra_matrix.getNumColumn(); j++) {
							f.setElem(i, j, 0);
						}
					}
					std::vector<int> ptr(extra_matrix.getNumLine());
					std::vector<int> d(extra_matrix.getNumLine());*/
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

/*bool Graph::checkContractionMinmality() {
	bool res = false;
	int maxFlow = 0;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numLine - 1; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1) {
				Graph new_graph(numNode - 2, k - 1);
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

				res = new_graph.algorithmEven();
				if (res == true) {
					return false;
				}
			}
		}
	}
	return true;
}*/

bool Graph::perebor_vertex_fun(int &num, std::vector<int> &index) {
	bool res = false, res_comp = true;
	int maxFlow = 0;
	num = k - 1;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numNode; i++) {
		if (!index.empty()) {
			res_comp = true;
			for (int j = 0; j < index.size(); j++) { //совпадает ли этот индекс с предыдущими?
				if (i == index[j]) {
					res_comp = false; //если да, то неверно
					break;
				}
			}
			if (res_comp == true) { //ели не совпадает, то делай дальше
				if (index.size() + 1 == num) { //если удалено нужное количество вершин, то начинай проверку
					Graph new_graph(numNode - num, k);
					Matrix new_matrix(numNode - num, numNode - num);
					int pos_i = 0, pos_j = 0;
					for (int k = 0; k < numNode; k++) { //построение нового графа
						for (int j = 0; j < index.size(); j++) {
							if (k == index[j] || k != i) {
								res_comp = false;
								break;
							}
						}
						if (res_comp) {
							for (int l = 0; l < numNode; l++) {
								for (int j = 0; j < index.size(); j++) {
									if (l == index[j] || l != i) {
										res_comp = false;
										break;
									}
								}
								if (res_comp) {
									new_matrix.setElem(pos_i, pos_j, matrix.getElem(k, l));
									pos_j++;
								}
							}
							pos_j = 0;
							pos_i++;
						}
					}
					new_graph.set_Matrix(new_matrix);

					res = new_graph.algorithmEven();
					if (res == true) {
						return false;
					}
				}
			}
		}
		if (res_comp == true) {
			index.push_back(i);
			perebor_vertex_fun(num, index);
			int id = index.size() - 1;
			index.erase(index.begin() + id);
		}
	}
	return true;
}

bool Graph::perebor_vertex_check() {
	bool res = false;
	int maxFlow = 0;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numNode; i++) {
		for (int j = 0; j < numNode; j++) {
			if (i != j) {
				for (int t = 0; t < numNode; t++) {
					if (i != t && j != t) {
						for (int w = 0; w < numNode; w++) {
							if (i != w && j != w && t != w) {
								Graph new_graph(numNode - 4, k);
								Matrix new_matrix(numNode - 4, numNode - 4);
								int pos_i = 0, pos_j = 0;
								for (int k = 0; k < numNode; k++) {
									if (k != i && k != j && k != t && k != w) {
										for (int l = 0; l < numNode; l++) {
											if (l != i && l != j && l != t && l != w) {
												new_matrix.setElem(pos_i, pos_j, matrix.getElem(k, l));
												pos_j++;
											}
										}
										pos_j = 0;
										pos_i++;
									}
								}
								new_graph.set_Matrix(new_matrix);

								res = new_graph.algorithmEven();
								if (res == true) {
									return false;
								}
							}
						}
					}
				}
			}
		}
	}
	return true;
}

bool Graph::connect(int s) {
	std::queue<int> q; // очередь с вершинами, которые мы рассматриваем на данном этапе
	q.push(s);
	std::vector<bool> used(numNode); // логический массив, указывающий, посещена ли вершина
	used[s] = true;
	while (!q.empty()) // пока мы не обойдем все вершины, которые можно достигнуть из данной
	{
		int v = q.front();
		q.pop(); // достаем из головы очереди одну вершину
		for (size_t i = 0; i < numLine; ++i) //просмотрим все ребра, исходящие из данной вершины
		{
			if (matrix.getElem(v, i) > 0) {
				int to = matrix.getElem(v, i);
				if (!used[to]) // если текущая вершина еще не была посещена
				{
					used[to] = true; //отмечаем, что мы ее посетили
					q.push(to); // помещаем в очередь
				}
			}
		}
	}
	std::vector<bool>::iterator it;
	it = find(used.begin(), used.end(), false); // проверяем, остались ли еще непосещенные вершины
	if (it == used.end())
		return true; // если все вершины посещены, то граф связный
	else
		return false;
}

bool Graph::perebor_vertex() {
	bool res = false;
	int maxFlow = 0;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numNode; i++) {
		for (int j = 0; j < numNode; j++) {
			if (i != j) {
				for (int t = 0; t < numNode; t++) {
					if (i != t && j != t) {
						Graph new_graph(numNode - 3, k);
						Matrix new_matrix(numNode - 3, numNode - 3);
						int pos_i = 0, pos_j = 0;
						for (int k = 0; k < numNode; k++) {
							if (k != i && k != j && k != t) {
								for (int l = 0; l < numNode; l++) {
									if (l != i && l != j && l != t) {
										new_matrix.setElem(pos_i, pos_j, matrix.getElem(k, l));
										pos_j++;
									}
								}
								pos_j = 0;
								pos_i++;
							}
						}
						new_graph.set_Matrix(new_matrix);

						for (int yu = 0; yu < new_graph.get_numLine(); yu++) {
							res = new_graph.connect(yu);
							if (res == true)
								std::cout << "res in connect = true" << std::endl;
							else
								std::cout << "res in connect = false" << std::endl;
							//res = new_graph.algorithmEven();
							//if (res == true) {
							if (res == false) {
								return false;
							}
						}
					}
				}
			}
		}
	}
	return true;
}

/*bool Graph::perebor_vertex() {
	bool res = false;
	int maxFlow = 0;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numNode; i++) {
		for (int j = 0; j < numNode; j++) {
			if (i != j) {
				
						Graph new_graph(numNode - 2, k);
						Matrix new_matrix(numNode - 2, numNode - 2);
						int pos_i = 0, pos_j = 0;
						for (int k = 0; k < numNode; k++) {
							if (k != i && k != j) {
								for (int l = 0; l < numNode; l++) {
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

						res = new_graph.algorithmEven();
						if (res == true) {
							return false;
						}

			}
		}
	}
	return true;
}*/

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

/*bool Graph::checkContractionMinmality() {
	bool res = false;
	int maxFlow = 0;
	if (numLine <= 0 || numColumn <= 0 || numNode <= 0) {
		return false;
	}
	for (int i = 0; i < numLine; i++) {
		for (int j = i + 1; j < numColumn; j++) {
			if (matrix.getElem(i, j) == 1) {
				Graph new_graph(numNode - 1, k);
				Matrix new_matrix(numNode - 1, numNode - 1);
				int I = 0, J = 0;
				int pos_i = 0, pos_j = 0;
				if (i < j) {
					I = i;
					J = j;
				}
				else {
					I = j;
					J = i;
				}
				for (int k = 0; k < numLine; k++) {
					if (k != I && k != J) {
						if (k < J) {
							pos_i = k;
						}
						else {
							pos_i = k - 1;
						}
						for (int l = 0; l < numColumn; l++) {
							if (l != I && l != J) {
								if (l < J) {
									pos_j = l;
								}
								else {
									pos_j = l - 1;
								}
								new_matrix.setElem(pos_i, pos_j, matrix.getElem(k, l));
							}
						}
						if (matrix.getElem(k, I) == 1 || matrix.getElem(k, J) == 1) {
							new_matrix.setElem(pos_i, I, 1);
						}
					}
				}
				for (int l = 0; l < numColumn; l++) {
					if (l != I && l != J) {
						if (l < J) {
							pos_j = l;
						}
						else {
							pos_j = l - 1;
						}
						if (matrix.getElem(I, l) == 1 || matrix.getElem(J, l) == 1)
							new_matrix.setElem(I, pos_j, 1);
					}
				}
				new_graph.set_Matrix(new_matrix);

				res = new_graph.algorithmEven();
				if (res == true) {
					return false;
				}
			}
		}
	}
	return true;
}*/

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