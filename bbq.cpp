#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <limits.h>
#include <unordered_map>
#include <time.h>
#include <fstream>
#include <stack>
//#include "queue_.h"
#include <pthread.h>
#include <algorithm>
#include <queue>
#include <string>
#include <sstream>
#include <ctime>

using namespace std;
static const int INF8 = 1000000000; // For unreachable pairs
size_t thread_count;
int k;

struct TreeRoot {
	    int father;
	    std::vector<uint32_t> vspt;
	    //std::vector<uint32_t> vmap; //the spt of v1  -> v[1]
	    //std::vector<std::vector<uint32_t> > d; // distance table
	    std::vector<uint32_t> child;
	    unordered_map<string, int>  d_map;
    };

    struct TreeNode {
	    int father;
	    std::vector<uint32_t> child;
	    std::vector<uint32_t> vspt;
	    //std::vector<std::vector<uint32_t> > d;
	    unordered_map<string, int>  d_map;
    };

    struct Tree {
	    TreeRoot root_;
	    std::vector<TreeNode> node_;
	    uint32_t height;
    };

    Tree tree_;
    int num_v_;
    std::vector<int> node_in_bag;

    double bbq_time[7];

    bool root_finish_flag = false;
    bool getdis(std::vector<std::unordered_map<int, uint32_t> > &index_, int v, int w) {
	    /*if (root_finish_flag == true)
	    	if ((tree_.root_.vmap[v] != uint32_t(-1)) && (tree_.root_.vmap[w] != uint32_t(-1))) {
	    		return tree_.root_.d[tree_.root_.vmap[v]][tree_.root_.vmap[w]];
	    	}
	    */
	    if (v == w) return 0;
	    std::unordered_map<int, uint32_t>::iterator iter = index_[v].find(w);
	    return (iter != index_[v].end()) ? iter->second : INT_MAX;
    }

    string i_to_string(int ival1, int ival2) {
	    stringstream ss1, ss2;
	    ss1 << ival1;
	    string s1 = ss1.str();

	    //string s2;
	    ss2 << ival2;
	    string s2 = ss2.str();

	    string s3 = s1 + "-" + s2;
	    return s3;
    }


    bool get_time(double *time) {
	    for (int i = 0; i < 7; i++) time[i] = bbq_time[i];
	    return true;
    }

    double get_time_pre() { return (bbq_time[0] + bbq_time[1] + bbq_time[2] + bbq_time[3] + bbq_time[4]); }

    vector<vector<int> > res;
    vector<int> sort_res;
    vector<int> spt_sort_res;
    pthread_mutex_t busy_lock = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t finished_lock = PTHREAD_MUTEX_INITIALIZER;

	vector<bool> busy;
	vector<bool> finished;

    void * b_to_t_apsp(void * rank);

int main(int argc, char** argv)
{
    if (argc != 3) {
        cout << "don't forget the value of k and thread " << endl;
    }

    k = atoi(argv[1]);
	//k = 2;
    thread_count = atoi(argv[2]);

	for (int i = 0; i < 7; i++) 
		bbq_time[i] = 0.0;

    
	ifstream ifs;
	ifs.open("rn-data.txt");
	ifs >> num_v_;
	ifs.close();
	int &V = num_v_;
	//tree_.root_.vmap.resize(V);
	//for (int i = 0; i < V; i++) tree_.root_.vmap[i] = -1;
	std::vector<int> deg(V);
	ifs.open("rn-data-d.txt");
	for (int v = 0; v < V; v++) ifs >> deg[v];
	ifs.close();
	std::vector<std::vector<int> > adj_v_(V);//adj: v
	{
		int tmp;
		ifs.open("rn-data-t-t.txt");

		for (int v = 0; v < V; v++) {
			for (int i = 0; i < deg[v]; i++) {
				ifs >> tmp;
				adj_v_[v].push_back(tmp - 1);
			}
		}
		ifs.close();
	}

	std::vector<std::unordered_map<int, uint32_t> > dis(V);//map -> key: v  value: distance
	vector<vector<int> > adj_d_(V); //adj:d

	{
		int tmp;
		ifs.open("rn-data-t-val.txt");
		for (int v = 0; v < V; ++v) {
			for (int i = 0; i < deg[v]; ++i) {
				ifs >> tmp;
				adj_d_[v].push_back(tmp);
			}
		}
		ifs.close();
		for (int v = 0; v < V; ++v) {
			for (int i = 0; i < deg[v]; ++i) {
				dis[v].insert(make_pair(adj_v_[v][i], adj_d_[v][i]));
			}
		}
	}
	


	double time1 = 0.0 - clock();
	std::vector<std::vector<int> > stack_;
	{
		//remove up to k
		//bool flag[V];//flag[v] = true: v has been removed.
		std::vector<bool> flag(V);
		for (int i = 0; i < V; i++) flag[i] = false;
		for (int degree = 1; degree <= k; degree++) {
			if (V - stack_.size() <= k + 1) break;
			for (int v = 0; v < V; v++) {
				if (V - stack_.size() <= k + 1) break;
				if ((flag[v] == false) && (adj_v_[v].size() <= degree)) {
					flag[v] = true;
					{
						std::vector<int> tmp;
						stack_.push_back(tmp);
					}
					stack_[stack_.size() - 1].push_back(v);
					for (size_t i = 0; i < adj_v_[v].size(); i++) {
						int v1 = adj_v_[v][i];
						stack_[stack_.size() - 1].push_back(v1);
						//remove v from v1
					std:
						vector<int>::iterator it = adj_v_[v1].begin();
						//vector<uint32_t>::iterator it1 = adj_d_[v1].begin();
						while (it != adj_v_[v1].end()) {
							if (*it == v) break;
							it++;
						}
						if (it != adj_v_[v1].end()) adj_v_[v1].erase(it);
						//if (it1 != adj_d_[v1].end()) adj_d_[v1].erase(it1);
						//adj_v_[v1].erase(remove(adj_v_[v1].begin(),adj_v_[v1].end(), v), adj_v_[v1].end());
					}	
				}
			}
		}
	/*	//root
		for (int v = 0; v < V; v++) {
			if (flag[v] == false) {
				tree_.root_.vmap[v] = (int)tree_.root_.vspt.size();
				tree_.root_.vspt.push_back(v);
			}
		}
	*/
		for (int v = 0; v < V; v++) {
			if (flag[v] == false) {
				tree_.root_.vspt.push_back(v);
			}
		}

		// 初始化根节点
		for (size_t i = 0; i < tree_.root_.vspt.size(); ++i) {
			int v1 = tree_.root_.vspt[i];
			for (size_t j = 0; j < tree_.root_.vspt.size(); ++j) {
				int v2 = tree_.root_.vspt[j];
				//tree_.root_.d[i].push_back(getdis, v1, v2);
				string str_v1_v2 = i_to_string(v1, v2);
				tree_.root_.d_map.insert(make_pair(str_v1_v2.c_str(), getdis(dis, v1, v2)));
			}
		}
	}
	bbq_time[0] = time1 + clock();
	std::ofstream out_prepro;
	out_prepro.open("out_of_preprocessing_2.txt");
	out_prepro << "graph -> stack " << bbq_time[0] / 1000000 << endl;
	//cout << "graph -> stack " << (clock() + time_preprocessing) / 1000 << endl;
	out_prepro << "root size = " << tree_.root_.vspt.size() << endl;
	//time_out << "graph -> stack " << (clock() + time_preprocessing) / 1000 << endl;
	//time_out << "root size = " << tree_.root_.vspt.size() << endl;
	//
	// local APSP: top(root)->bottom
	//
	//root APSP
	
	//int time_out_flag = 5000;
	time1 = 0.0 - clock();


/*
	tree_.root_.d.resize(tree_.root_.vspt.size());
	std::vector<uint32_t> dist_to(V);
	
	
	for (size_t si = 0; si < tree_.root_.vspt.size(); si++) {
		for (int v = 0; v < V; v++) dist_to[v] = INF8;
		int s = tree_.root_.vspt[si];
		dist_to[s] = 0;
		index_min_pq pq(V);
		pq.insert(s, 0);
		while (!pq.isempty()) {
			int v = pq.delmin();
			int dist_s_v = dist_to[v];
			for (size_t i = 0; i < adj_v_[v].size(); i++) {
				int w = adj_v_[v][i];
				int ew = adj_d_[v][i];
				if (dist_to[w] > dist_s_v + ew) {
					dist_to[w] = dist_s_v + ew;
					if (pq.contains(w)) pq.update(w, dist_s_v + ew);
					else pq.insert(w, dist_s_v + ew);
				}
			}
		}

		for (size_t i = 0; i < tree_.root_.vspt.size(); i++) {
			int vi = tree_.root_.vspt[i];
			//if (getdis(dis, s, vi) == INT_MAX) dis[s].insert(std::make_pair(vi, dist_to[vi]));
			//else update_value(dis, s, vi, dist_to[vi]);
			tree_.root_.d[si].push_back(dist_to[vi]);
		}
		//pq.~index_min_pq();
		pq.del_index();
		pq.~index_min_pq(); 
	}
*/

	
	bbq_time[1] = time1 + clock();
	out_prepro << "root dij " << bbq_time[1] / 1000 << endl;
	//cout << "root dij " << (clock() + time1) / 1000 << endl;
	//time_out << "root dij " << (clock() + time1) / 1000 << endl;
	//root_finish_flag = true;


	time1 = 0.0 - clock();
/*
	for (size_t bag = stack_.size() - 1; bag < stack_.size(); bag--) {
		int v = stack_[bag][0];
		for (size_t i1 = 1; i1 < stack_[bag].size() - 1; i1++) {
			int v1 = stack_[bag][i1];
			for (size_t i2 = i1 + 1; i2 < stack_[bag].size(); i2++) {
				int v2 = stack_[bag][i2];
				int dis_v1_v2 = getdis(dis, v1, v2);
				int tmp = getdis(dis, v, v1) + getdis(dis, v, v2);
				if (dis_v1_v2 < tmp) {//edge(v1, v2) may change the distance
					for (size_t i3 = 0; i3 < stack_[bag].size(); i3++)
					{
						if ((i3 == i1) || (i3 == i2)) continue;
						int v3 = stack_[bag][i3];
						int dis_v3_v1 = getdis(dis, v3, v1);
						int dis_v3_v2 = getdis(dis, v3, v2);
						if (dis_v3_v1 < dis_v3_v2) {//may update dis_v3_v2
							if (dis_v3_v2 > dis_v3_v1 + dis_v1_v2)
								update_value(dis, v3, v2, dis_v3_v1 + dis_v1_v2);
						}
						else {//may update dis_v3_v1
							if (dis_v3_v1 > dis_v3_v2 + dis_v1_v2)
								update_value(dis, v3, v1, dis_v3_v2 + dis_v1_v2);
						}
					}
                }
			}
		}
	}
*/

	bbq_time[2] = time1 + clock();
	out_prepro << "bag dij " <<  bbq_time[2]/ 1000000 << endl;
	//cout << "bag dij " << (clock() + time1) / 1000 << endl;
	//time_out << "bag dij " << (clock() + time1) / 1000 << endl;
	//
	// stack -> tree
	//
	time1 = 0.0 - clock();

	
	{
		tree_.node_.resize(stack_.size()+1);
		node_in_bag.resize(num_v_);
		for (size_t i = 0; i < num_v_; i++)
			node_in_bag[i] = -2;//have not add to the tree
		for (size_t i = 0; i < tree_.root_.vspt.size(); i++)
			node_in_bag[tree_.root_.vspt[i]] = -1;//node is in root
		std::vector<int> rm_order(num_v_);//record the order of node removing
		for (size_t i = 0; i < tree_.root_.vspt.size(); i++)
			rm_order[tree_.root_.vspt[i]] = (int)stack_.size();
		for (size_t i = 0; i < stack_.size(); i++)
			rm_order[stack_[i][0]] = (int)i;
		
		tree_.root_.father = -1;
		tree_.node_[0].father = tree_.root_.father;
		tree_.node_[0].child.assign(tree_.root_.child.begin(), tree_.root_.child.end());
		tree_.node_[0].vspt.assign(tree_.root_.vspt.begin(), tree_.root_.vspt.end());
		//tree_.node_[0].d(tree_.root_.d);
		//tree_.node_[0].d_map(tree_.root_.d_map);
		for (size_t i = 0; i < tree_.root_.vspt.size(); ++i) {
			int v1 = tree_.root_.vspt[i];
			for (size_t j = 0; j < tree_.root_.vspt.size(); ++j) {
				int v2 = tree_.root_.vspt[j];
				//tree_.root_.d[i].push_back(getdis, v1, v2);
				string str_v1_v2 = i_to_string(v1, v2);
				tree_.node_[0].d_map.insert(make_pair(str_v1_v2.c_str(), getdis(dis, v1, v2)));
			}
		}

		

		for (size_t bag = stack_.size() - 1; bag < stack_.size(); bag--) {
			//get the first removing node
			int v_first = stack_[bag][1];
			for (size_t i = 2; i < stack_[bag].size(); i++) {
				int v_tmp = stack_[bag][i];
				if (rm_order[v_tmp] < rm_order[v_first])
					v_first = v_tmp;
			}

			size_t nodespt = stack_.size() - bag - 1;
			node_in_bag[stack_[bag][0]] = (int)nodespt;
			if (rm_order[v_first] == stack_.size()) {
				tree_.node_[nodespt+1].father = 0;
				tree_.node_[0].child.push_back((uint32_t)nodespt+1);
			
			}
			else {
				tree_.node_[nodespt+1].father = node_in_bag[v_first] + 1;


				tree_.node_[node_in_bag[v_first] + 1].child.push_back((uint32_t)nodespt+1);



			}
			for (size_t i = 0; i < stack_[bag].size(); i++)
				tree_.node_[nodespt+1].vspt.push_back(stack_[bag][i]);
			//tree_.node_[nodespt+1].d.resize(stack_[bag].size());

			
			for (size_t i = 0; i < stack_[bag].size(); i++) {
				int v1 = stack_[bag][i];
				for (size_t j = 0; j < stack_[bag].size(); j++) {
					int v2 = stack_[bag][j];
					//tree_.node_[nodespt+1].d[i].push_back(getdis(dis, v1, v2));
					string str_v1_v2 = i_to_string(v1, v2);
					tree_.node_[nodespt+1].d_map.insert(make_pair(str_v1_v2.c_str(), getdis(dis, v1, v2)));
				}
			}
		}
	}

	
    bbq_time[3] = time1 + clock();
	out_prepro << "stack -> tree " << bbq_time[3] / 1000000 << endl;
	out_prepro.close();

    std::ofstream out1;
	out1.open("tree.txt");
	out1 << "the riginal root:" << endl;
	for (size_t i = 0; i < tree_.root_.vspt.size(); i++) out1 << tree_.root_.vspt[i] << " ";
	out1 << endl;
	out1 << "the new tree:" << endl;
	for (size_t node = 0; node < tree_.node_.size(); node++) {
		for (size_t v = 0; v < tree_.node_[node].vspt.size(); v++) out1 << tree_.node_[node].vspt[v] << " ";
		out1 << endl;
	}

    //得到树的自底向上的层次遍历的结果

    //vector<vector<int> > res;
	queue<int> q;
	q.push(0);
	
	while (!q.empty()) {
		//vector<TreeNode> one_level;
		vector<int> one_level;

		size_t size_ = q.size();
		for (int i = 0; i < size_; ++i) {
			int tem_num = q.front();
			q.pop();
			one_level.push_back(tem_num);
			for (int j = 0; j < tree_.node_[tem_num].child.size(); ++j) {
				q.push(tree_.node_[tem_num].child[j]);
			}
		}
		res.push_back(one_level);
	}
	reverse(res.begin(), res.end());//从底向上的遍历结果
	
	//vector<int> sort_res;
	
	for (int i = 0; i < res.size(); ++i) {
		for (int j = 0; j < res[i].size(); ++j) {
			sort_res.push_back(res[i][j]);
		}
	}

    //将树的节点号与存储它的vector的下标对应
	size_t len_sort_res = sort_res.size();
	for (size_t i = 0; i < len_sort_res; ++i) {
		spt_sort_res.push_back(0);
	}

	for (size_t i = 0; i < len_sort_res; ++i) {
		busy.push_back(false);
	}

	for (size_t i = 0; i < len_sort_res; ++i) {
		finished.push_back(false);
	}

	//vector<int> spt_sort_res(len_sort_res, 0);
	for (size_t i = 0; i < sort_res.size(); ++i) {
		spt_sort_res[sort_res[i]] = i;
	}
	
    
    size_t thread;
	pthread_t * thread_handles;
	//thread_handles = malloc(thread_count * sizeof(pthread_t));
	thread_handles = new pthread_t[thread_count];
	//double start_c = clock();
	//double end_c;
	for ( thread = 0; thread < thread_count; ++thread ){
		pthread_create(&thread_handles[thread], NULL, b_to_t_apsp, (void*) thread );
	}

	//cout << "Hello from the main thread" << endl;
	//end_c = clock();
	//double time_last = end_c - start_c;
	//cout << "the time from buttom to top: " << time_last / 1000000 << endl;

	for (thread = 0; thread < thread_count; ++thread) {
		pthread_join(thread_handles[thread], NULL);
	}
	//end_c = clock();
	//double time_last = end_c - start_c;
	//cout << "the time from buttom to top: " << (end_c - start_c) / 1000000 << endl;

    delete [] thread_handles;
    return 0;
}

//线程函数
void * b_to_t_apsp(void * rank) {
		
	long myrank = (long) rank;
	time_t thread_sum_time = 0;
	//double start_c;
	time_t t_start = time(NULL);
	
	for (int i = 0; i < sort_res.size() ; ++i) {
		pthread_mutex_lock(&busy_lock);
		if (busy[i] == false)
			busy[i] = true;
		//pthread_mutex_unlock(&busy_lock);
		else {
			pthread_mutex_unlock(&busy_lock);
			continue;
		}
		//double start_c = clock();
		pthread_mutex_unlock(&busy_lock);
		time_t t_start = time(NULL);
		//double start_c = clock();
		//vector<int> sort_node;
		//sort_node.assign(tree_.node_[sort_res[i]].vspt.begin(), tree_.node_[sort_res[i]].vspt.end());
		//int v0 = sort_node[0];

			
		//int v0 = tree_.node_[sort_res[i]].vspt[0];
		for (uint32_t j = 0; j <  tree_.node_[sort_res[i]].vspt.size() - 1; ++j) {
			int v1 =  tree_.node_[sort_res[i]].vspt[j];
			for (uint32_t k = j + 1; k < tree_.node_[sort_res[i]].vspt.size(); ++k) {
				int v2 =  tree_.node_[sort_res[i]].vspt[k];
				string thr_str_1_2 = i_to_string(v1, v2);
				//string thr_str_0_1 = i_to_string(v1, v0);
				//string thr_str_0_2 = i_to_string(v0, v2);
				int dis_v1_v2 = tree_.node_[sort_res[i]].d_map[thr_str_1_2.c_str()];
				//先将父亲节点中的每条边都与他的每个孩子节点的中的边比较，更新父亲节点
				uint32_t child_num = tree_.node_[sort_res[i]].child.size();
				//int count = 0;
				for (uint32_t m = 0; m < child_num; ++m) {
					uint32_t  child_number = tree_.node_[sort_res[i]].child[m];
					uint32_t spt_child_m = spt_sort_res[child_number];
						
					unordered_map<string, int>::iterator map_it = tree_.node_[child_number].d_map.find(thr_str_1_2.c_str());
					if (map_it != tree_.node_[child_number].d_map.end()) {
						pthread_mutex_lock(&finished_lock);
						while (!finished[spt_child_m]); 
						pthread_mutex_unlock(&finished_lock);
							if ( map_it->second < dis_v1_v2 )
								tree_.node_[sort_res[i]].d_map[thr_str_1_2.c_str()] = map_it->second;
							
						
					}
					
					
				}
				
			}
		}
		int v0 = tree_.node_[sort_res[i]].vspt[0];
		for (uint32_t j = 1; j <  tree_.node_[sort_res[i]].vspt.size() - 1; ++j) {
			int v1 =  tree_.node_[sort_res[i]].vspt[j];
			for (uint32_t k = j + 1; k < tree_.node_[sort_res[i]].vspt.size(); ++k) {
				int v2 =  tree_.node_[sort_res[i]].vspt[k];
				string thr_str_1_2 = i_to_string(v1, v2);
				string thr_str_0_1 = i_to_string(v1, v0);
				string thr_str_0_2 = i_to_string(v0, v2);
				int dis_v1_v2 = tree_.node_[sort_res[i]].d_map[thr_str_1_2.c_str()];
				int dis_tmp = tree_.node_[sort_res[i]].d_map[thr_str_0_1.c_str()] + tree_.node_[sort_res[i]].d_map[thr_str_0_2.c_str()];
				if (dis_v1_v2 > dis_tmp) {
					tree_.node_[sort_res[i]].d_map[thr_str_1_2.c_str()] = dis_tmp;
				}		
				
			}

		
			pthread_mutex_lock(&finished_lock);
			finished[i] = true;
			pthread_mutex_unlock(&finished_lock);
		
		}		

	}
	//double end_c = clock();
	time_t t_end = time(NULL);
    thread_sum_time += difftime(t_end, t_start);
	//thread_sum_time += (end_c-start_c) / CLOCKS_PER_SEC;		
	cout << " thread " << myrank  << " finished" << " ......" << " time: " 
		<< thread_sum_time << " s" << endl;
	return NULL;
}
	
	
	







