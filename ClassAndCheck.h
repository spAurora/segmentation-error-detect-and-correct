//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////
#ifndef _CLASSANDCHECK_H_
#define _CLASSANDCHECK_H_

using namespace std;
using namespace cv;

class CRegion
{
public:
	int id;   //区域编号
	vector<int> pixelLocation;   //区域的各个像素位置

	int pixelNum;   //像素数（面积）

	CRegion()  //无参构造
	{
		id = -1;
		pixelNum = 0;
	}
protected:
private:
};

class BTNode
{
public:
	BTNode* left;
	BTNode* right;
	BTNode* father;

	bool hasMerged;
	bool highLevelNode;
	bool anchored;

	int ID;
	int area;
	int level;
	vector<int> pixelLocation;
	double avgB;
	double avgG;
	double avgR;

	vector<int> basicNodeID;

	double spectualStandDeviation_B;
	double spectualStandDeviation_G;
	double spectualStandDeviation_R;
	int boundingBoxPerimeter;
	int perimeter;

	BTNode()
	{
		left = NULL;
		right = NULL;
		father = NULL;
		hasMerged = false;
		anchored = false;
		highLevelNode = true;
		ID = -1;
		area = -1;
		avgB = -1;
		avgG = -1;
		avgR = -1;
		
		spectualStandDeviation_B = -1;
		spectualStandDeviation_G = -1;
		spectualStandDeviation_R = -1;
		boundingBoxPerimeter = -1;
		perimeter = -1;
	}
	//重载+运算符
	//仅统计面积以及合并像素
	BTNode operator+(const BTNode & b)  
	{
		BTNode bTNode;
		bTNode.area = this->area + b.area;
		for (int i = 0; i<this->pixelLocation.size(); i++)
			bTNode.pixelLocation.push_back(this->pixelLocation[i]);
		for (int i = 0; i<b.pixelLocation.size(); i++)
			bTNode.pixelLocation.push_back(b.pixelLocation[i]);
		return bTNode;
	}
};

class GraphNode
{
public:
	int ID; 
	GraphNode()
	{
		ID = -1;
	}
	GraphNode(int mID)
	{
		ID = mID;
	}
};

class ArrayHeadGraphNode  //头结点数组
{
public:
	vector<GraphNode> pGraphNodeList;		//邻接拓扑点
protected:
private:
};

class Edge
{
public:
	int key;
	int n_1;
	int n_2;
	double length;
};

void CheckRegionNum(int regionNum)
{
	printf("-----------\n");
	printf("check regionNum:%d\n", regionNum);
}

void CheckRegionSet(CRegion* cRegion)
{
	printf("-----------\n");
	printf("check region set:\nregion_0 pixelNum:%d\n", cRegion[0].pixelLocation.size());
}

void CheckGplot(ArrayHeadGraphNode* head)
{
	printf("-----------\n");
	printf("check gplot:\n");
	vector<GraphNode>::iterator it;
	for (it = head[0].pGraphNodeList.begin(); it!= head[0].pGraphNodeList.end(); it++)
		printf("%d -> ", it->ID);
	printf("\n");
	for (it = head[1].pGraphNodeList.begin(); it!= head[1].pGraphNodeList.end(); it++)
		printf("%d -> ", it->ID);
	printf("\n");
	for (it = head[2].pGraphNodeList.begin(); it!= head[2].pGraphNodeList.end(); it++)
		printf("%d -> ", it->ID);
	printf("\n");
}

void CheckBTNode(BTNode* bTNode)
{
	printf("-----------\n");
	printf("check bTNode:\n");
	for (int i = 0; i<3; i++)
	{
		printf("AVG BGR:%lf %lf %lf\n", bTNode[i].avgB, bTNode[i].avgG, bTNode[i].avgR);
	}
}

void CheckScaleSetsInitialization(vector<Edge> & edge, ArrayHeadGraphNode* head, BTNode* bTNode)
{
	printf("-----------\n");
	printf("check initialize:\n");
	printf("edge vector size: %d\n", edge.size());
	printf("length: %lf->%lf->%lf\n", edge[0].length, edge[1].length, edge[2].length);
	printf("head[0]:");
	vector<GraphNode>::iterator it;
	for (it = head[0].pGraphNodeList.begin(); it!= head[0].pGraphNodeList.end(); it++)
		printf("%d -> ", it->ID);
	printf("\n");
	printf("bTnode[0] ");
	printf("AVG BGR:%lf %lf %lf\n",  bTNode[0].avgB, bTNode[0].avgG, bTNode[0].avgR);
	printf("spectual deviation: %lf %lf %lf\n", bTNode[0].spectualStandDeviation_B, bTNode[0].spectualStandDeviation_G, bTNode[0].spectualStandDeviation_R);
	printf("Bounding Box Perimeter: %d\n", bTNode[0].boundingBoxPerimeter);
	printf("Perimeter:%d\n", bTNode[0].perimeter);
}

void CheckScaleSetsBuildingProcess(int N, int i)
{
	int TenPercent = (N-1)/10;
	if (i == 0)
	{
		printf("-----------\n");
		printf("start scale-sets building...\n");
	}
		if (i % TenPercent == 0 && i != 0)
	{
		printf("merge progress:%d%%...\n", i/TenPercent*10);
	}
}

void CheckScaleSets(vector<Edge> & edge, ArrayHeadGraphNode* head, BTNode* bTNode, int regionNum)
{
	printf("-----------\n");
	printf("check scale sets:\n");
	printf("now edge size: %d\n", edge.size());
	printf("bTNode final node area:%d\n", bTNode[2*regionNum-2].area);
	printf("-----------\n");
	for (int i = regionNum; i<regionNum+3; i++)
	{
		printf("merge region check: %d\n", i-regionNum);
		printf("ID：%d\nchildren ID:%d %d\n",bTNode[i].ID, bTNode[i].left->ID, bTNode[i].right->ID);
		printf("adjacency list：\n");
		for (int j = 0; j<head[i].pGraphNodeList.size(); j++)
		{
			printf("%d->", head[i].pGraphNodeList[j].ID);
		}
		printf("\n");
	}
}

int CheckNewRegionNum(int* labels, int width, int height)
{
	printf("-----------\n");
	printf("check new regionNum:\n");
	int max = 0;
	vector<int> num;
	for (int i = 0; i<height; i++)
		for (int j = 0; j<width; j++)
		{
			bool exist = false;
			if (labels[i*width+j] > max)
				max = labels[i*width+j];
			for (int k = 0; k < num.size(); k++)
			{
				if (labels[i*width+j] == num[k])
				{
					exist = true;
					break;
				}
			}
			if (exist == false)
			{
				num.push_back(labels[i*width + j]);
			}
		}
	printf("new regionNum: %d\n", max);
	printf("new regionNum + 1: %d\n", num.size());
	return max;
}

void CheckVKMIFirstFinal(double VK_first, double VK_final, double MI_first, double MI_final)
{
	printf("-----------\n");
	printf("check VK MI first&final: ");
	printf("%lf, %lf, %lf, %lf\n", VK_first, VK_final, MI_first, MI_final);
}

void CheckConflictRegionNum(int conflictNum)
{
	printf("-----------\n");
	printf("check conflict region num: %d\n", conflictNum);
}

void CheckMergeErrorDetentionStart()
{
	printf("-----------\n");
	printf("merge error deteneion has started...\n");
}

void CheckWhichBeAnchored(bool* mainIsBetter, int size)
{
	printf("-----------\n");
	printf("Which Be Anchored List:\n");
	int count = 0;
	for (int i = 0; i<size; i++)
	{
		printf("%d, ", mainIsBetter[i]);
		if (mainIsBetter[i] == false)
			count++;
	}
	printf("\n");
	printf("Total conflict region pairs number: %d\n", size);
	printf("Anchor region num:%d\n", count);

}

void CheckLoopNum(int t)
{
	printf("\n\n***********\n");
	printf("NOW LOOP NUM: %d\n", t);
	printf("***********\n");

}

void CheckEdge(vector<Edge> & edge)
{
	printf("-----------\n");
	printf("check edge size: %d\n", edge.size());
	for (int i = 0; i<edge.size(); i++)
	{
		printf("%d, %d\n", edge[i].n_1, edge[i].n_2);
	}
	printf("\n");
}

#endif