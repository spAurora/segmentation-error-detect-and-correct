//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <list>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <forward_list>
#include <cmath>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "highgui.h"
#include "cv.h"


#include "ScaleSets.h"
#include "Output.h"
#include "Evaluation.h"
#include "MergeErrorDetention.h"
#include "test.h"

using namespace std;
using namespace cv;

#pragma comment(linker, "/STACK:102400000,102400000")    //防止栈溢出

int main()
{
	//Test_1();
	int numSuperpixels = 3000;  //初始分割的区域数
	double compactness = 10;   //初始分割的紧凑度
	int Q = 6;  //欠分割惩罚，值越大，全局尺度越小
	int bestLocationMode = 2;  //0是每次迭代重新计算最优尺度位置,1是利用Q固定最优尺度，2是人工输入最优尺度
	double relative_location = 0.88;   //人工最优尺度模式下的相对最优尺度

	char load_srimg[] = "./data/dqh-0.6m-4.bmp";
	char output_path[] = "./output/";
	char evalution_path[] = "./output/evalution.txt";

	int width, height, regionNum;
	Mat srimg;
	//读取原始影像
	LoadSrimg(srimg, load_srimg);

	width = srimg.cols;
	height = srimg.rows;	
	
	vector<BTNode> newAnchorMap;
	vector<int> anchoredRegionID;
	vector<BTNode> bTNode_LevelOne;
	int t = 0;

	FILE *fp;
	if((fp = fopen(evalution_path, "w+")) == NULL)
	{
		printf("打开evalution记录文件失败\n");
		system("pause");
		exit(-1);
	}
	
	/********************************主循环start*******************************/
	int regionNum_lock_c1 = -1, regionNum_lock_c2 = -1;
	do 
	{
		CheckLoopNum(++t);
		int* labels = new int[height*width];
		if (newAnchorMap.size() == 0)
		{
			//初始分割
			SlicSegmentation(width, height, compactness, numSuperpixels, regionNum, srimg, labels);
		}
		else
		{
			//储存被锚定的区域ID
			anchoredRegionID.clear();
			for (int j = 0; j<newAnchorMap.size(); j++)
				if (newAnchorMap[j].anchored == true)
					anchoredRegionID.push_back(j);
			//记录区域数
			regionNum = newAnchorMap.size();
			//填充连通图
			for (int j = 0; j<newAnchorMap.size(); j++)
				for (int k =0; k<newAnchorMap[j].pixelLocation.size(); k++)
				{
					labels[newAnchorMap[j].pixelLocation[k]] = j;
				}
		}
		CheckRegionNum(regionNum);
		CRegion* cRegion = new CRegion[regionNum];
		ArrayHeadGraphNode* head = new ArrayHeadGraphNode[regionNum];
		//初始化区域集合
		CreateRegionSet(labels, srimg, cRegion, regionNum, width, height);
		CheckRegionSet(cRegion);
		//初始化拓扑图
		CreateToplogicalGraph(labels, head, regionNum, width, height);
		CheckGplot(head);
		BTNode* bTNode_basic = new BTNode[regionNum];
		//初始化二元分割树底层节点
		CreateBasicBTNodeArray(bTNode_basic, srimg, cRegion, regionNum, anchoredRegionID);
		CheckBTNode(bTNode_basic);
		//如果是第一次计算则保存初始分割结果
		if (t == 1)
		{
			for (int i = 0 ;i<regionNum; i++)
				bTNode_LevelOne.push_back(bTNode_basic[i]);
		}
		//输出初始分割结果
		if (t == 1)
			OutputSegmentResult(labels, srimg, output_path, "FirstOverSegmentResult.bmp");
		delete[] cRegion;
		delete[] labels;

		//构建尺度集
		ArrayHeadGraphNode* head_c1 = new ArrayHeadGraphNode[2*regionNum - 1];
		BTNode* bTNode_c1 = new BTNode[2*regionNum - 1];
		CreateScaleSetsModel(head, head_c1, bTNode_basic, bTNode_c1, regionNum, CalcuteCriterion_HeterogeneityChange, srimg);

		ArrayHeadGraphNode* head_c2 = new ArrayHeadGraphNode[2*regionNum - 1];
		BTNode* bTNode_c2 = new BTNode[2*regionNum - 1];
		CreateScaleSetsModel(head, head_c2, bTNode_basic, bTNode_c2, regionNum, CalcuteCriterion_MaxSpectralDifference, srimg);



		//计算最优尺度并输出分割结果
		int regionNum_LevelOne = bTNode_LevelOne.size();
		BTNode* bTNode_LevelOneDArray = new BTNode[regionNum_LevelOne];
		for (int i = 0; i<regionNum_LevelOne; i++)
			bTNode_LevelOneDArray[i] = bTNode_LevelOne[i];
		int best_location_c1, best_location_c2;

		if (bestLocationMode == 1)
		{
			if (regionNum_lock_c1 == -1)
				best_location_c1 = CalculateBestScaleLocation(bTNode_c1, regionNum, srimg, Q, bTNode_LevelOneDArray, regionNum_LevelOne);
			else
				best_location_c1 = regionNum - regionNum_lock_c1 - 1;
		}
		else if (bestLocationMode == 0)
		{
			best_location_c1 = CalculateBestScaleLocation(bTNode_c1, regionNum, srimg, Q, bTNode_LevelOneDArray, regionNum_LevelOne);
		}
		else if (bestLocationMode == 2)
		{
			if (regionNum_lock_c1 == -1)
				best_location_c1 = regionNum * relative_location;
			else
				best_location_c1 = regionNum - regionNum_lock_c1 - 1;
		}
		printf("best location c1:%d\n", best_location_c1);
		int* outLabels_avgSpec = new int[height*width];
		SetAllNodeValueByAbsoluteLocation(outLabels_avgSpec, bTNode_c1, regionNum, best_location_c1);
		char merge_result_main[256];
		sprintf(merge_result_main, "MergeResultMain_%d.bmp", t);
		OutputSegmentResult(outLabels_avgSpec, srimg, output_path, merge_result_main);
		if (regionNum_lock_c1 == -1)
		{
			regionNum_lock_c1 = CheckNewRegionNum(outLabels_avgSpec, width, height);
		}
		else
			CheckNewRegionNum(outLabels_avgSpec, width, height);
		delete[] outLabels_avgSpec;

		if (bestLocationMode == 1)
		{
			if (regionNum_lock_c2 == -1)
				best_location_c2 = CalculateBestScaleLocation(bTNode_c2, regionNum, srimg, Q, bTNode_LevelOneDArray, regionNum_LevelOne);
			else
				best_location_c2 = regionNum - regionNum_lock_c2 - 1;
		}
		else if (bestLocationMode == 0)
		{
			best_location_c2 = CalculateBestScaleLocation(bTNode_c2, regionNum, srimg, Q, bTNode_LevelOneDArray, regionNum_LevelOne);
		}
		else if (bestLocationMode == 2)
		{
			if (regionNum_lock_c2 == -1)
				best_location_c2 = regionNum * relative_location;
			else
				best_location_c2 = regionNum - regionNum_lock_c2 - 1;
		}
		printf("best location c2:%d\n", best_location_c2);
		int* outLabels_maxSpec = new int[height*width];
		SetAllNodeValueByAbsoluteLocation(outLabels_maxSpec, bTNode_c2, regionNum, best_location_c2);
		char merge_result_ref[256];
		sprintf(merge_result_ref, "MergeResultRef_%d.bmp", t);
		OutputSegmentResult(outLabels_maxSpec, srimg, output_path, merge_result_ref);
		if (regionNum_lock_c2 == -1)
		{
			regionNum_lock_c2 = CheckNewRegionNum(outLabels_maxSpec, width, height);
		}
		else
			CheckNewRegionNum(outLabels_maxSpec, width, height);
		delete[] outLabels_maxSpec;

		//保存最优尺度下的区域集合
		int regionNum_main = regionNum - best_location_c1;
		int regionNum_ref = regionNum - best_location_c2;
		BTNode* bTNode_main = new BTNode[regionNum_main];
		BTNode* bTNode_ref = new BTNode[regionNum_ref];
		int location_main = 0, location_ref = 0;
		for (int i = 0; i<regionNum + best_location_c1; i++)
			if (bTNode_c1[i].highLevelNode == true)
			{
				bTNode_main[location_main] = bTNode_c1[i];
				location_main++;
			}
		for (int i = 0; i<regionNum + best_location_c2; i++)
			if (bTNode_c2[i].highLevelNode == true)
			{
				bTNode_ref[location_ref] = bTNode_c2[i];
				location_ref++;
			}
		if (location_main != regionNum_main || location_ref != regionNum_ref)
			{
				printf("location main or ref wrong!");
					system("pause");
			}

		//输出量化评估结果
		double e_main, e_ref;
		//e_main = Evalution_E(bTNode_main, regionNum_main, srimg);
		//e_ref = Evalution_E(bTNode_ref, regionNum_ref, srimg);
		e_main = Evalution_Z(bTNode_main, regionNum_main, srimg);
		e_ref = Evalution_Z(bTNode_ref, regionNum_ref, srimg);
		printf("\nooooooooooooooooooooooooooooooooooooooooooooo/n%d %d %.3lf %.3lf\n\nooooooooooooooooooooooooooooooooooooooooooooo\n",regionNum_main, regionNum_ref, e_main, e_ref);
		fprintf(fp, "%d %d %.3lf %.3lf\n",regionNum_main, regionNum_ref, e_main, e_ref);


		//清空newAnchorMap准备填充
		newAnchorMap.clear();
		int existAnchor;
		existAnchor = MergeErrorDetention_Main(bTNode_basic, regionNum, bTNode_main, bTNode_ref, regionNum_main, regionNum_ref, srimg, newAnchorMap);
		//若无冲突跳出循环
		if (existAnchor == 0)
			break;
		int *outLabels_newAnchorMap = new int[width*height];
		SetLabelsValueByBTNodeList(newAnchorMap, outLabels_newAnchorMap);

		//输出锚定图
		char anchormap[256], anchormap_highlight[256], anchormap_fill[256];
		sprintf(anchormap, "AnchorMap_%d.bmp", t);
		sprintf(anchormap_highlight, "AnchorMapHighLight_%d.bmp", t);
		sprintf(anchormap_fill, "AnchorMapFill_%d.bmp", t);
		OutputSegmentResult(outLabels_newAnchorMap, srimg, output_path, anchormap);
		OutputAnchorMap(outLabels_newAnchorMap, srimg, output_path, anchormap_highlight, newAnchorMap);
		OutputAnchorMap_fill(outLabels_newAnchorMap, srimg, output_path, anchormap_fill, newAnchorMap);

		delete[] bTNode_c1;
		delete[] bTNode_c2;
		delete[] head;
		delete[] head_c1;
		delete[] head_c2;
		delete[] bTNode_LevelOneDArray;
		delete[] bTNode_main;
		delete[] bTNode_ref;
		delete[] bTNode_basic;
	} while (1);
	/********************************主循环end***********************************/

	fclose(fp);
	system("pause");
	return 0;
}