//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////
#ifndef _SCALESETS_H_
#define _SCALESETS_H_

#include "Slic.h"
#include "ClassAndCheck.h"
#include "Criterion.h"

#define min(a, b) ((a)<(b)?(a):(b))
#define max(a,b,c) (((((a)>(b)) ? (a):(b))>(c)) ? (((a)>(b)) ? (a):(b)):(c))
#define MINC 9999999999999

using namespace std;          
using namespace cv;

void LoadSrimg(Mat & srimg, char load_img[])
{	
	/*
	*读取原始影像
	*/
	srimg = imread(load_img, 1);
	if (srimg.empty())
	{
		printf("Can not open Image\n");
		system("pause");
		exit(0);
	}
}

void SlicSegmentation(int width, int height, int compactness, int numSuperpixels, int & regionNum, Mat srimg, int* labels)
{
	/*
	*SLIC超像素分割
	*生成底层过分割区域
	*/
	int step, numseeds;
	int sz = width*height;
	int* rin = new int[height*width];
	int* gin = new int[height*width];
	int* bin = new int[height*width];
	double* lvec = new double[height*width];
	double* avec = new double[height*width];
	double* bvec = new double[height*width];
	int* klabels = new int[height*width];
	int* seedIndices = new int[height*width];
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
		{
			bin[i*width + j] = srimg.data[(i*width + j)*3];
			gin[i*width + j] = srimg.data[(i*width + j)*3 + 1];
			rin[i*width + j] = srimg.data[(i*width + j)*3 + 2];
		}
		//RGB空间转换为LAB空间
		rgbtolab(rin,gin,bin,sz,lvec,avec,bvec);
		//寻找种子点
		step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5;
		getLABXYSeeds(step,width,height,seedIndices,&numseeds);
		double* kseedsx = new double[numseeds];
		double* kseedsy = new double[numseeds];
		double* kseedsl = new double[numseeds];
		double* kseedsa = new double[numseeds];
		double* kseedsb = new double[numseeds];
		for(int k = 0; k < numseeds; k++)
		{
			kseedsx[k] = seedIndices[k]%width;
			kseedsy[k] = seedIndices[k]/width;
			kseedsl[k] = lvec[seedIndices[k]];
			kseedsa[k] = avec[seedIndices[k]];
			kseedsb[k] = bvec[seedIndices[k]];
		}
		//计算超像素
		PerformSuperpixelSLIC(lvec, avec, bvec, kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,width,height,numseeds,klabels,step,compactness);
		//生成连通图
		EnforceSuperpixelConnectivity(klabels, width, height, numSuperpixels, labels, &regionNum);
}

void CreateRegionSet(int* labels, Mat &srimg, CRegion* cRegion, int regionNum, int width, int height)
{
	/*
	*创建区域集合
	*/
	for (int i = 0; i<regionNum; i++)
		cRegion[i].id = i;
	for (int i = 0;i<height;i++)
		for (int j = 0;j<width;j++)
		{
			cRegion[labels[i*width + j]].pixelLocation.push_back(i*width+j);
			cRegion[labels[i*width + j]].pixelNum++;
		}
}


bool Cmp_head(GraphNode first, GraphNode second) 
{
	/*
	*sort函数的第三个参数
	*将head的邻接节点根据ID从小到大排序
	*/
	return first.ID < second.ID;
}


void CreateToplogicalGraph(int* labels, ArrayHeadGraphNode* head, int regionNum, int width, int height)
{
	/*
	*创建区域的拓扑图
	*邻接结点排序
	*/
	for (int i = 0; i<height-1; i++)
		for (int j = 0; j<width-1; j++)
		{
			if (labels[i*width + j] != labels[i*width + j + 1])    //横向
			{
				vector<GraphNode>::iterator it;
				int check = 0; //0为不存在
				for (it = head[labels[i*width + j]].pGraphNodeList.begin(); it != head[labels[i*width + j]].pGraphNodeList.end(); it++)
					if (it->ID == labels[i*width+j+1])
					{
						check = 1;
						break;
					}
				if (check == 0)
				{
					head[labels[i*width + j]].pGraphNodeList.push_back(labels[i*width + j + 1]);
					head[labels[i*width + j + 1]].pGraphNodeList.push_back(labels[i*width + j]);
				}
			}

			if (labels[i*width + j] != labels[(i+1)*width+j])  //纵向
			{
				vector<GraphNode>::iterator it;
				int check = 0; //0为不存在
				for (it = head[labels[i*width + j]].pGraphNodeList.begin(); it != head[labels[i*width + j]].pGraphNodeList.end(); it++)
					if (it->ID == labels[(i+1)*width+j])
					{
						check = 1;
						break;
					}
				if (check == 0)
				{
					head[labels[i*width + j]].pGraphNodeList.push_back(labels[(i+1)*width + j]);
					head[labels[(i+1)*width + j]].pGraphNodeList.push_back(labels[i*width + j]);
				}
			}
		}
	for(int i = 0; i<regionNum; i++)
		sort(head[i].pGraphNodeList.begin(), head[i].pGraphNodeList.end(), Cmp_head);
}

void CreateBasicBTNodeArray(BTNode* btNode, Mat srimg, CRegion* cRegion, int regionNum, vector<int> & anchoredRegionID)
{
	/*
	*初始化二元分割树底层节点
	*/
	for(int i = 0;  i<regionNum;  i++)
	{
		btNode[i].ID = cRegion[i].id;
		btNode[i].area = cRegion[i].pixelNum;
		btNode[i].pixelLocation.assign(cRegion[i].pixelLocation.begin(), cRegion[i].pixelLocation.end());
		btNode[i].basicNodeID.push_back(i);
		btNode[i].avgB = 0;
		btNode[i].avgG = 0;
		btNode[i].avgR = 0;
		for(int j = 0; j<btNode[i].pixelLocation.size(); j++)
		{
			btNode[i].avgB += srimg.data[btNode[i].pixelLocation[j]*3];
			btNode[i].avgG += srimg.data[btNode[i].pixelLocation[j]*3+1];
			btNode[i].avgR += srimg.data[btNode[i].pixelLocation[j]*3+2];
		}
		btNode[i].avgB /= btNode[i].area;
		btNode[i].avgG /= btNode[i].area;
		btNode[i].avgR /= btNode[i].area;
	}
	//标记已经被锚定的区域
	for (int i = 0; i<anchoredRegionID.size(); i++)
		btNode[anchoredRegionID[i]].anchored = true;
}


bool IsBothMinEdge(ArrayHeadGraphNode* head, int sID, int eID, BTNode* bTNode, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*判断是否为双向最小边（闭合环）
	*/
	int min_regionID = -1;
	double min_c = MINC;
	for (int i = 0; i<head[sID].pGraphNodeList.size(); i++)
	{
		double temp_c = (*F)(&bTNode[sID], &bTNode[head[sID].pGraphNodeList[i].ID], srimg);
		if (temp_c < min_c)
		{
			min_c = temp_c;
			min_regionID = head[sID].pGraphNodeList[i].ID;
		}
	}
	if (min_regionID == eID)
		return true;
	else
		return false;
}

bool HasConstructedClosedLoop(int id, vector<Edge> & edge)
{
	/*
	*判断是否参与构建闭合环
	*/
	for(int i = 0; i<edge.size(); i++)
		if(edge[i].n_2 == id)
			return true;
	return false;
}

bool Cmp_edge(Edge first, Edge second) 
{
	/*
	*sort函数的第三个参数
	*将edge的根据length小到大排
	*/
	return first.length < second.length;
}

void CreateEdgeSequenceList(vector<Edge> & edge, ArrayHeadGraphNode* head, int regionNum, BTNode* bTNode, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*初始化边集合
	*/
	for (int i =0; i<regionNum; i++)
	{
		//检查是否为锚定区域
		if (bTNode[i].anchored == true)
			continue;
		//检查是否已参与构建闭合环
		if (HasConstructedClosedLoop(i, edge))
			continue;
		//寻找最小边
		double min_c = MINC;
		int min_regionID = -1;
		for (int j = 0;  j<head[i].pGraphNodeList.size();  j++)
		{
			//首先检查是否被锚定
			if (bTNode[head[i].pGraphNodeList[j].ID].anchored == true)
				continue;

			double temp_c = (*F)(&bTNode[i], &bTNode[head[i].pGraphNodeList[j].ID], srimg);
			if (temp_c < min_c)
			{
				min_c = temp_c;
				min_regionID = head[i].pGraphNodeList[j].ID;
			}
		}
		//检查是否为双向最小边
		bool isBothMinEdge = false;
		if (min_regionID != -1)
		{
			isBothMinEdge = IsBothMinEdge(head, min_regionID, i, bTNode, F, srimg);
		}
		else
			isBothMinEdge = false;
		//是双向最小边将该边加入边顺序表
		if (isBothMinEdge == true)
		{
			Edge temp_edge;
			temp_edge.length = min_c;
			temp_edge.n_1 = i;
			temp_edge.n_2 = min_regionID;
			edge.push_back(temp_edge);
		}
	}
	//根据边长从小到大排序
	sort(edge.begin(), edge.end(), Cmp_edge);
}


void InitializeScaleSetsModel(vector<Edge> & edge, ArrayHeadGraphNode* head_basic, ArrayHeadGraphNode* head, BTNode* bTNode_basic, BTNode* bTNode, int regionNum, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*尺度集模型初始化
	*包括双向最短边、拓扑图、二元分割树节点共3部分
	*/
	CreateEdgeSequenceList(edge, head_basic, regionNum, bTNode_basic, F, srimg);
	for (int i = 0; i<regionNum; i++)
	{
		bTNode[i] = bTNode_basic[i];
		head[i] = head_basic[i];
	}
}

void UpdateBinaryPartitionTree(BTNode* bTNode, int location, int childID_1, int childID_2, int regionNum, Mat & srimg)
{
	/*
	*更新二元分割树
	*将合并后的新节点加入二元分割树中
	*属性更新应保证满足异质性准则的计算
	*/
	location += regionNum;
	bTNode[location].ID = location;
	bTNode[location].left = &bTNode[childID_1];
	bTNode[location].right = &bTNode[childID_2];
	bTNode[location].left->father = &bTNode[location];
	bTNode[location].right->father = &bTNode[location];
	bTNode[location].area = bTNode[childID_1].area + bTNode[childID_2].area;
	bTNode[location].avgB = (bTNode[childID_1].area*bTNode[childID_1].avgB + bTNode[childID_2].area*bTNode[childID_2].avgB) / bTNode[location].area;
	bTNode[location].avgG = (bTNode[childID_1].area*bTNode[childID_1].avgG + bTNode[childID_2].area*bTNode[childID_2].avgG) / bTNode[location].area;
	bTNode[location].avgR = (bTNode[childID_1].area*bTNode[childID_1].avgR + bTNode[childID_2].area*bTNode[childID_2].avgR) / bTNode[location].area;
	for (int i = 0; i<bTNode[childID_1].pixelLocation.size(); i++)
		bTNode[location].pixelLocation.push_back(bTNode[childID_1].pixelLocation[i]);
	for (int i = 0; i<bTNode[childID_2].pixelLocation.size(); i++)
		bTNode[location].pixelLocation.push_back(bTNode[childID_2].pixelLocation[i]);
	for (int i = 0; i<bTNode[childID_1].basicNodeID.size(); i++)
		bTNode[location].basicNodeID.push_back(bTNode[childID_1].basicNodeID[i]);
	for (int i = 0; i<bTNode[childID_2].basicNodeID.size(); i++)
		bTNode[location].basicNodeID.push_back(bTNode[childID_2].basicNodeID[i]);
}

void UpdateAdjacentRegion(ArrayHeadGraphNode* head, int location, int inID, int delId_1, int delId_2)
{
	/*
	更新邻接区域
	*/
	vector<GraphNode>::iterator it = head[inID].pGraphNodeList.begin();
	//删除参与合并的两区域
	while(it != head[inID].pGraphNodeList.end())
	{
		if (it->ID == delId_1 || it->ID == delId_2)
			//该写法不会发生迭代器越界
			it = head[inID].pGraphNodeList.erase(it);
		else
			it++;
	}
	//插入合并生成的新区域
	head[inID].pGraphNodeList.push_back(GraphNode(location));
}

void UpdateTopologicalGraph(ArrayHeadGraphNode* head, int location, int childID_1, int childID_2, int regionNum)
{
	/*
	*更新邻接拓扑图
	*参与合并的两个区域的邻接区域取并集，同时更新所有邻接区域
	*保持邻接区域ID的有序性可以极大提升计算速度
	*/
	location += regionNum;
	vector<GraphNode>::iterator it_1 = head[childID_1].pGraphNodeList.begin();
	vector<GraphNode>::iterator it_2 = head[childID_2].pGraphNodeList.begin();
	while(it_1 != head[childID_1].pGraphNodeList.end() && it_2 != head[childID_2].pGraphNodeList.end())
	{
		//参与合并的两个区域不进入并集
		if (it_1->ID == childID_1 || it_1->ID == childID_2)
		{
			it_1++;
			continue;
		}
		if (it_2->ID == childID_1 || it_2->ID == childID_2)
		{
			it_2++;
			continue;
		}
		if (it_1->ID < it_2->ID)
		{
			head[location].pGraphNodeList.push_back(*it_1);
			it_1++;
		}
		else if (it_2->ID < it_1->ID)
		{
			head[location].pGraphNodeList.push_back(*it_2);
			it_2++;
		}
		else
		{
			head[location].pGraphNodeList.push_back(*it_1);
			it_1++;
			it_2++;
		}
	}
	while(it_1 != head[childID_1].pGraphNodeList.end())
	{
		if (it_1->ID != childID_1 && it_1->ID != childID_2)
		{
			head[location].pGraphNodeList.push_back(*it_1);
			it_1++;
		}
		else
			it_1++;
	}
	while(it_2 != head[childID_2].pGraphNodeList.end())
	{
		if (it_2->ID != childID_1 && it_2->ID != childID_2)
		{
			head[location].pGraphNodeList.push_back(*it_2);
			it_2++;
		}
		else
			it_2++;
	}
	//更新邻接区域
	vector<GraphNode>::iterator it;
	for (it = head[location].pGraphNodeList.begin(); it!= head[location].pGraphNodeList.end(); it++)
		UpdateAdjacentRegion(head, location, it->ID, childID_1, childID_2);
}

void UpdateEdge(vector<Edge> & edge, int location, int regionNum, ArrayHeadGraphNode* head, BTNode* bTNode, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*更新边（闭合环）集合
	*/
	//location += regionNum;
	////删除已经参与合并的最短边
	//edge.erase(edge.begin());
	////更新所有邻接区域的闭合环信息
	//for (int i = 0; i<head[location].pGraphNodeList.size();i++)
	//{
	//	int sID = head[location].pGraphNodeList[i].ID;
	//	double min_c = MINC;
	//	int min_regionID = -1;
	//	if (bTNode[sID].anchored == false)
	//	{
	//		for (int j = 0; j<head[sID].pGraphNodeList.size(); j++)
	//		{
	//			if (bTNode[head[sID].pGraphNodeList[j].ID].anchored == false)
	//			{
	//				double temp_c = (*F)(&bTNode[sID], &bTNode[head[sID].pGraphNodeList[j].ID], srimg);
	//				if (temp_c < min_c)
	//				{
	//					min_c = temp_c;
	//					min_regionID = head[sID].pGraphNodeList[j].ID;
	//				}
	//			}
	//		}
	//		bool isBothMinEdge = false;
	//		//if (min_regionID != -1)
	//		//{
	//			isBothMinEdge = IsBothMinEdge(head, min_regionID, i, bTNode, F, srimg);
	//		//}
	//		//else
	//			//isBothMinEdge = false;
	//		if (isBothMinEdge == false)
	//		{
	//			//若无闭合环，只需要在edge中删除包含相关区域的边即可
	//			vector<Edge>::iterator it;
	//			//注意下面的写法，如果不这么写迭代器可能会越界
	//			for (it = edge.begin(); it != edge.end();)
	//				if (it->n_1 == sID || it->n_2 == sID)
	//					it = edge.erase(it);
	//				else
	//					it++;
	//		}
	//		else
	//		{
	//			//若存在闭合环,需要检查重复
	//			//若有重复跳过，若无重复需要删除相关再插入新闭合环
	//			bool rEpeat = false;
	//			for (int k = 0; k<edge.size(); k++)
	//				if ((edge[k].n_1 == sID && edge[k].n_2 == min_regionID) || (edge[k].n_2 == sID && edge[k].n_1 == min_regionID))
	//					rEpeat = true;
	//			if (rEpeat == false)
	//			{
	//				//删除包含相关区域的边
	//				vector<Edge>::iterator it = edge.begin();
	//				while (it != edge.end())
	//				{
	//					if (it->n_1 == sID || it->n_2 == sID || it->n_1 == min_regionID || it->n_2 == min_regionID)
	//						it = edge.erase(it);
	//					else
	//						it++;
	//				}
	//				//插入新闭合环
	//				Edge temp_edge;
	//				temp_edge.length = min_c;
	//				temp_edge.n_1 = sID;
	//				temp_edge.n_2 = min_regionID;
	//				vector<Edge>::iterator itt;
	//				if (edge.empty())
	//					edge.push_back(temp_edge);
	//				else
	//					for (itt = edge.begin(); itt != edge.end(); itt++)
	//					{
	//						if (itt == edge.end())
	//						{
	//							edge.insert(itt, temp_edge);
	//							break;
	//						}
	//						if (itt->length >= min_c)
	//						{
	//							edge.insert(itt, temp_edge);
	//							break;
	//						}
	//					}
	//			}
	//		}
	//	}
	//	
	//}
	location += regionNum;
	//删除已经参与合并的最短边
	edge.erase(edge.begin());
	//更新所有邻接区域的闭合环信息
	for (int i = 0; i<head[location].pGraphNodeList.size();i++)
	{
		int sID = head[location].pGraphNodeList[i].ID;
		if (bTNode[sID].anchored == true)
			continue;
		double min_c = MINC;
		int min_regionID = -1;
		for (int j = 0; j<head[sID].pGraphNodeList.size(); j++)
		{
			if (bTNode[head[sID].pGraphNodeList[j].ID].anchored == true)
				continue;
			double temp_c = (*F)(&bTNode[sID], &bTNode[head[sID].pGraphNodeList[j].ID], srimg);
			if (temp_c < min_c)
			{
				min_c = temp_c;
				min_regionID = head[sID].pGraphNodeList[j].ID;
			}
		}
		//bool isBothMinEdge = IsBothMinEdge(head, min_regionID, sID, bTNode, F, srimg);
		bool isBothMinEdge = false;
		if (min_regionID != -1)
			isBothMinEdge = IsBothMinEdge(head, min_regionID, sID, bTNode, F, srimg);
		else
			isBothMinEdge = false;
		if (isBothMinEdge == false)
		{
			//若无闭合环，只需要在edge中删除包含相关区域的边即可
			vector<Edge>::iterator it;
			//注意下面的写法，如果不这么写迭代器可能会越界
			for (it = edge.begin(); it != edge.end();)
				if (it->n_1 == sID || it->n_2 == sID)
					it = edge.erase(it);
				else
					it++;
		}
		else
		{
			//若存在闭合环,需要检查重复
			//若有重复跳过，若无重复需要删除相关再插入新闭合环
			bool rEpeat = false;
			for (int k = 0; k<edge.size(); k++)
				if ((edge[k].n_1 == sID && edge[k].n_2 == min_regionID) || (edge[k].n_2 == sID && edge[k].n_1 == min_regionID))
					rEpeat = true;
			if (rEpeat == false)
			{
				//删除包含相关区域的边
				vector<Edge>::iterator it = edge.begin();
				while (it != edge.end())
				{
					if (it->n_1 == sID || it->n_2 == sID || it->n_1 == min_regionID || it->n_2 == min_regionID)
						it = edge.erase(it);
					else
						it++;
				}
				//插入新闭合环
				Edge temp_edge;
				temp_edge.length = min_c;
				temp_edge.n_1 = sID;
				temp_edge.n_2 = min_regionID;
				vector<Edge>::iterator itt;
				if (edge.empty())
					edge.push_back(temp_edge);
				else
					for (itt = edge.begin(); itt != edge.end(); itt++)
					{
						if (itt == edge.end())
						{
							edge.insert(itt, temp_edge);
							break;
						}
						if (itt->length >= min_c)
						{
							edge.insert(itt, temp_edge);
							break;
						}
					}
			}
		}
	}
}

void CreateEdgeSequenceList_new(vector<Edge> & edge, ArrayHeadGraphNode* head, int empty_location, BTNode* bTNode, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*初始化边集合(合并的第二阶段)
	*/
	for (int i = 0; i < empty_location; i++)
	{
		//检查是否为锚定区域
		if (bTNode[i].anchored == true)
			continue;
		//检查是否已经参与合并
		if (bTNode[i].hasMerged == true)
			continue;
		//寻找最小边
		double min_c = MINC;
		int min_regionID = -1;
		for (int j = 0;  j<head[i].pGraphNodeList.size();  j++)
		{
			//首先检查是否被锚定
			if (bTNode[head[i].pGraphNodeList[j].ID].anchored == true)
				continue;
			//然后检查是否已参与合并
			if (bTNode[head[i].pGraphNodeList[j].ID].hasMerged == true)
				continue;

			double temp_c = (*F)(&bTNode[i], &bTNode[head[i].pGraphNodeList[j].ID], srimg);
			if (temp_c < min_c)
			{
				min_c = temp_c;
				min_regionID = head[i].pGraphNodeList[j].ID;
			}
		}
		//不需要检查是否是双向最小边
		//但必须保证单向最短边存在
		if (min_regionID != -1)
		{
			Edge temp_edge;
			temp_edge.length = min_c;
			temp_edge.n_1 = i;
			temp_edge.n_2 = min_regionID;
			edge.push_back(temp_edge);
		}
	}
	//根据边长从小到大排序
	sort(edge.begin(), edge.end(), Cmp_edge);
}

void UpdateEdge_new(vector<Edge> & edge, int location, int regionNum, ArrayHeadGraphNode* head, BTNode* bTNode, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*更新单向最短边集合（第二阶段）
	*/
	location += regionNum;
	//删除参与合并区域的最短边（共2条）
	int n_1 = edge[0].n_1;
	int n_2 = edge[0].n_2;
	vector<Edge>::iterator it_0;
	for (it_0 = edge.begin(); it_0 != edge.end();)
		if (it_0->n_1 == n_2)
			it_0 = edge.erase(it_0);
		else
			it_0++;
	edge.erase(edge.begin());
	//更新所有邻接区域的单向最短边信息
	for (int i = 0; i<head[location].pGraphNodeList.size();i++)
	{
		int sID = head[location].pGraphNodeList[i].ID;
		//跳过锚定的区域
		if (bTNode[sID].anchored == true)
			continue;
		//首先删除以该点为起点的原最短边
		vector<Edge>::iterator it;
		//注意下面的写法，如果不这么写迭代器可能会越界
		for (it = edge.begin(); it != edge.end();)
			if (it->n_1 == sID)
				it = edge.erase(it);
			else
				it++;
		double min_c = MINC;
		int min_regionID = -1;
		for (int j = 0; j<head[sID].pGraphNodeList.size(); j++)
		{
			//跳过锚定的区域
			if (bTNode[head[sID].pGraphNodeList[j].ID].anchored == true)
				continue;
			double temp_c = (*F)(&bTNode[sID], &bTNode[head[sID].pGraphNodeList[j].ID], srimg);
			if (temp_c < min_c)
			{
				min_c = temp_c;
				min_regionID = head[sID].pGraphNodeList[j].ID;
			}
		}
		if (min_regionID != -1)
		{
			Edge temp_edge;
			temp_edge.length = min_c;
			temp_edge.n_1 = sID;
			temp_edge.n_2 = min_regionID;
			vector<Edge>::iterator itt;
			if (edge.empty())
				edge.push_back(temp_edge);
			else
				for (itt = edge.begin(); itt != edge.end(); itt++)
				{
					if (itt == edge.end())
					{
						edge.insert(itt, temp_edge);
						break;
					}
					if (itt->length >= min_c)
					{
						edge.insert(itt, temp_edge);
						break;
					}
				}
		}
	}
}

void CreateScaleSetsModel(ArrayHeadGraphNode* head_basic, ArrayHeadGraphNode* head, BTNode* bTNode_basic, BTNode* bTNode, int regionNum, double (*F)(BTNode*, BTNode*, Mat), Mat & srimg)
{
	/*
	*构建尺度集模型
	*/
	vector<Edge> edge;
	//初始化
	InitializeScaleSetsModel(edge, head_basic, head, bTNode_basic, bTNode, regionNum, F, srimg);
    CheckScaleSetsInitialization(edge, head, bTNode);
	//循环构建尺度集
	const int N = regionNum;
	bool emptyProcess = false;
	int emptyProcess_Location = 0;
	for (int i = 0; i<N-1; i++)
	{
		if (edge.empty())
		{
			printf("wrong!edge array empty!\n");
			return;
		}
		int n_1 = edge[0].n_1;
		int n_2 = edge[0].n_2;
		CheckScaleSetsBuildingProcess(N, i);
		//更新二元分割树
		UpdateBinaryPartitionTree(bTNode, i, n_1, n_2, regionNum, srimg);
		//更新拓扑图
		UpdateTopologicalGraph(head, i, n_1, n_2, regionNum);
		//更新双向最短边（闭合环）集合
		UpdateEdge(edge, i, regionNum, head, bTNode, F, srimg);
		//检查双向边集合是否为空
		if (edge.empty() && i != N-2)
		{
			printf("edge empty!\n");
			printf("edge empty location = %d\n", i);
			//system("pause");
			emptyProcess = true;
			emptyProcess_Location = i+1;
			break;
		}
	}
	//合并第二阶段，单向最短边合并
	if (emptyProcess == true)  
	{
		printf("merge part 2...\n");
		for (int i = regionNum; i<regionNum + emptyProcess_Location; i++)
		{
			bTNode[bTNode[i].left->ID].hasMerged = true;
			bTNode[bTNode[i].right->ID].hasMerged = true;
		}
		bool emptyProcess_new = false; 
		int newEmptyProcess_Location = -1;
		CreateEdgeSequenceList_new(edge, head, regionNum + emptyProcess_Location, bTNode, F, srimg);
		//CheckEdge(edge);
		for (int i = emptyProcess_Location; i<N-1; i++)
		{
			if (edge.empty())
			{
				printf("worng! edge array empty!\n");
				return;
			}
			int n_1 = edge[0].n_1;
			int n_2 = edge[0].n_2;
			//printf("mergeid: %d, %d, now edge size:%d\n", n_1, n_2, edge.size());

			CheckScaleSetsBuildingProcess(N, i);
			//更新二元分割树
			UpdateBinaryPartitionTree(bTNode, i, n_1, n_2, regionNum, srimg);
			//更新拓扑图
			UpdateTopologicalGraph(head, i, n_1, n_2, regionNum);
			//更新最短边集合
			UpdateEdge_new(edge, i, regionNum, head, bTNode, F, srimg);
			if (edge.empty() && i != N-2)
			{
				printf("edge empty(part 2)!\n");
				printf("edge empty location = %d\n", i);
				//system("pause");
				emptyProcess_new = true;
				newEmptyProcess_Location = i+1;
				break;
			}
			//int mergeID_1 = -2, mergeID_2 = -2;
			//bool existUnanchoredAdjacent = false;
			//for (int j = 0; j < regionNum + i; j++)
			//{
			//	if(bTNode[j].hasMerged == false && bTNode[j].anchored == false)
			//	{
			//		mergeID_1 = j;
			//		//检索未锚定的邻接区域
			//		for (int k = 0; k<head[j].pGraphNodeList.size(); k++)
			//		{
			//			if (bTNode[head[j].pGraphNodeList[k].ID].anchored == false)
			//			{
			//				mergeID_2 = head[j].pGraphNodeList[k].ID;
			//				existUnanchoredAdjacent = true;
			//				break;
			//			}
			//		}
			//		//目标区域存在未锚定的邻接区域
			//		if (mergeID_2 != -2)
			//		{
			//			bTNode[mergeID_1].hasMerged = true;
			//			bTNode[mergeID_2].hasMerged = true;
			//			break;
			//		}
			//		else
			//			//目标区域的邻接区域都被锚定了
			//			continue;
			//	}
			//}
			////存在未锚定的区域邻接
			//if (existUnanchoredAdjacent == true)
			//{
			//	UpdateBinaryPartitionTree(bTNode, i, mergeID_1, mergeID_2, regionNum, srimg);
			//	UpdateTopologicalGraph(head, i, mergeID_1, mergeID_2, regionNum);
			//}
			//else
			//{
			//	newEmptyProcess_Location = i;
			//	printf("no adjacent unanchored region pair!\n");
			//	printf("no adjacent unanchored region pair location = %d\n", i);
			//	break;
			//}
		}
		//CheckEdge(edge);
		//第三阶段，自由合并锚定区域
		if (emptyProcess_new == true)
		{
			printf("free merge(part 3)...!\n");
			for (int i = regionNum; i<regionNum + newEmptyProcess_Location; i++)
			{
				bTNode[bTNode[i].left->ID].hasMerged = true;
				bTNode[bTNode[i].right->ID].hasMerged = true;
			}
			for (int i = newEmptyProcess_Location; i<N-1; i++)
			{
				int mergeID_1, mergeID_2;
				for (int j = 0; j < regionNum + i; j++)
				{
					if(bTNode[j].hasMerged == false)
					{
						mergeID_1 = j;
						mergeID_2 = head[j].pGraphNodeList[0].ID;
						bTNode[mergeID_1].hasMerged = true;
						bTNode[mergeID_2].hasMerged = true;
						break;
					}
				}
				UpdateBinaryPartitionTree(bTNode, i, mergeID_1, mergeID_2, regionNum, srimg);
				UpdateTopologicalGraph(head, i, mergeID_1, mergeID_2, regionNum);
			}
		}
	}
	CheckScaleSets(edge, head, bTNode, regionNum);
}

#endif