//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////
#ifndef _MERGEERRORDETENTION_H_
#define _MERGEERRORDETENTION_H_

#include "ClassAndCheck.h"
#include "Evaluation.h"

using namespace std;
using namespace cv;

#define T 0.008

class CConflictPair
{
public:
	int location_main;
	int location_ref;

	vector<int> onlyMainNode;
	vector<int> onlyRefNode;
	vector<int> iNNode;

	CConflictPair()
	{
		location_main = -1;
		location_ref = -1;
	}
protected:
private:
};

bool IsBothConflict(BTNode bTNode_1, BTNode bTNode_2)
{
	/*
	�ж���������Ƿ��������ͻ
	*/
	for (int i = 0; i<bTNode_1.basicNodeID.size(); i++)
		for (int j = 0; j<bTNode_2.basicNodeID.size(); j++)
			if (bTNode_1.basicNodeID[i] == bTNode_2.basicNodeID[j])
				return true;
	return false;
}

int ConflictMode(BTNode & bTNode_main, BTNode & bTNode_ref)
{
	/*
	*�жϳ�ͻ����
	*/
	sort(bTNode_main.basicNodeID.begin(), bTNode_main.basicNodeID.end());
	sort(bTNode_ref.basicNodeID.begin(), bTNode_ref.basicNodeID.end());
	
	for (int i = 0; i<bTNode_main.basicNodeID.size(); i++)
	{
		bool cmp = false;
		for (int j = 0; j<bTNode_ref.basicNodeID.size(); j++)
		{
			if (bTNode_main.basicNodeID[i] == bTNode_ref.basicNodeID[j])
			{
				cmp = true;
				break;
			}
		}
		//���Ӽ�
		if (cmp == false)
			return 2;
	}
	if (bTNode_ref.basicNodeID.size() == bTNode_main.basicNodeID.size())
		//���
		return 0;   
	else
		//���Ӽ�
		return 1;   
}

void FindRegionConflict(vector<BTNode> & bTNode_basic, int regionNum_basic, BTNode* bTNode_main, BTNode* bTNode_ref, int regionNum_main, int regionNum_ref, vector<CConflictPair> & cConflictPair, Mat srimg)
{
	/*
	Ѱ�������ͻ
	*/
	for (int i = 0; i<regionNum_main; i++)
		for (int j = 0; j<regionNum_ref; j++)
			if (IsBothConflict(bTNode_main[i], bTNode_ref[j]) == true)
			{
				CConflictPair cConflictPair_temp;
				cConflictPair_temp.location_main = i;
				cConflictPair_temp.location_ref = j;
				//�б��ͻ����
				bool isMergeError = false;
				int conflictmode = -1;
				if (bTNode_main[i].basicNodeID.size() <= bTNode_ref[j].basicNodeID.size())
					conflictmode = ConflictMode(bTNode_main[i], bTNode_ref[j]);
				else
					conflictmode = ConflictMode(bTNode_ref[j], bTNode_main[i]);

                if (conflictmode == 2)
				{

					/*printf("conflict ID: %d, %d\n", i, j);
					printf("conflict mode: %d\n", conflictmode);
					printf("main: ");
					for (int k = 0; k<bTNode_main[i].basicNodeID.size(); k++)
						printf("%d->", bTNode_main[i].basicNodeID[k]);
					printf("\n");
					printf("ref: ");
					for (int k = 0; k<bTNode_ref[j].basicNodeID.size(); k++)
						printf("%d->", bTNode_ref[j].basicNodeID[k]);
					printf("\n");*/

					isMergeError = true;
					//��� ��Ӧ��ʼ�������������
					vector<int>::iterator it_1 = bTNode_main[i].basicNodeID.begin();
					vector<int>::iterator it_2 = bTNode_ref[j].basicNodeID.begin();
					while(it_1 != bTNode_main[i].basicNodeID.end() && it_2 != bTNode_ref[j].basicNodeID.end())
					{
						if (*it_1 < *it_2)
						{
							cConflictPair_temp.onlyMainNode.push_back(*it_1);
							it_1++;
							continue;
						}
						if (*it_1 > *it_2)
						{
							cConflictPair_temp.onlyRefNode.push_back(*it_2);
							it_2++;
							continue;
						}
						if (*it_1 == *it_2)
						{
							cConflictPair_temp.iNNode.push_back(*it_1);
							it_1++;
							it_2++;
							continue;
						}
					}
					//������
					while(it_1 != bTNode_main[i].basicNodeID.end())
					{
						cConflictPair_temp.onlyMainNode.push_back(*it_1);
						it_1++;
					}
					while(it_2 != bTNode_ref[j].basicNodeID.end())
					{
						cConflictPair_temp.onlyRefNode.push_back(*it_2);
						it_2++;
					}
					/*printf("only main node: ");
					for (int k = 0; k<cConflictPair_temp.onlyMainNode.size(); k++)
						printf("%d->", cConflictPair_temp.onlyMainNode[k]);
					printf("\n");
					printf("only ref node: ");
					for (int k = 0; k<cConflictPair_temp.onlyRefNode.size(); k++)
						printf("%d->", cConflictPair_temp.onlyRefNode[k]);
					printf("\n");
					printf("IN Node: ");
					for (int k = 0; k<cConflictPair_temp.iNNode.size(); k++)
						printf("%d->", cConflictPair_temp.iNNode[k]);
					printf("\n\n");*/
				}
                if (isMergeError == true)
					cConflictPair.push_back(cConflictPair_temp);
				//break;
			}
	//����main��������ref������ͻ��ɾ������ĳ�ͻ��
	printf("deleting redundant region pair...\n");
	vector<CConflictPair>::iterator it = cConflictPair.begin();
	while(it != cConflictPair.end())
	{
		vector<CConflictPair>::iterator it_later = it + 1;
		while(it_later != cConflictPair.end())
		{
			if (it_later->location_main == it->location_main)
				it_later = cConflictPair.erase(it_later);
			else
				it_later++;
		}
		it++;
	}
}

void BTNodeReInitialization(BTNode & bTNode, Mat srimg)
{
	/*
	*���³�ʼ��BTNode
	*/
	bTNode.area = bTNode.pixelLocation.size();
	bTNode.avgB = 0;
	bTNode.avgG = 0;
	bTNode.avgR = 0;
	for(int j = 0; j<bTNode.pixelLocation.size(); j++)
	{
		bTNode.avgB += srimg.data[bTNode.pixelLocation[j]*3];
		bTNode.avgG += srimg.data[bTNode.pixelLocation[j]*3+1];
		bTNode.avgR += srimg.data[bTNode.pixelLocation[j]*3+2];
	}
	bTNode.avgB /= bTNode.area;
	bTNode.avgG /= bTNode.area;
	bTNode.avgR /= bTNode.area;
}

bool MainIsBetterThanRef(vector<BTNode> & bTNode_AnchorMap, int regionNum_anchor, CConflictPair cConflictPair, Mat srimg)
{
	/*
	*�жϳ�ͻ����Main��Ref��һ������
	*/
	int mergeNum = cConflictPair.onlyMainNode.size() + cConflictPair.onlyRefNode.size() + cConflictPair.iNNode.size();
	int regionNum = regionNum_anchor - mergeNum + 2;
	for (int i = 0; i<regionNum_anchor; i++)
		bTNode_AnchorMap[i].hasMerged = false;

	for (int i = 0; i<cConflictPair.onlyMainNode.size(); i++)
		bTNode_AnchorMap[cConflictPair.onlyMainNode[i]].hasMerged = true;
	for (int i = 0; i<cConflictPair.iNNode.size(); i++)
		bTNode_AnchorMap[cConflictPair.iNNode[i]].hasMerged = true;
	for (int i = 0; i<cConflictPair.onlyRefNode.size(); i++)
		bTNode_AnchorMap[cConflictPair.onlyRefNode[i]].hasMerged = true;
	//���main
	BTNode* bTNode_AnchorMap_MainHold = new BTNode[regionNum];
	BTNode bTNode_main_temp = bTNode_AnchorMap[cConflictPair.onlyMainNode[0]];
	BTNode bTNode_ref_temp = bTNode_AnchorMap[cConflictPair.onlyRefNode[0]];
	for (int i = 1; i<cConflictPair.onlyMainNode.size(); i++)
		bTNode_main_temp = bTNode_main_temp + bTNode_AnchorMap[cConflictPair.onlyMainNode[i]];
	for (int i = 0; i<cConflictPair.iNNode.size(); i++)
		bTNode_main_temp = bTNode_main_temp + bTNode_AnchorMap[cConflictPair.iNNode[i]];
	for (int i = 1; i<cConflictPair.onlyRefNode.size(); i++)
		bTNode_ref_temp = bTNode_ref_temp + bTNode_AnchorMap[cConflictPair.onlyRefNode[i]];
	bTNode_AnchorMap_MainHold[regionNum - 1] = bTNode_main_temp;
	bTNode_AnchorMap_MainHold[regionNum - 1].ID = regionNum - 1;
	bTNode_AnchorMap_MainHold[regionNum - 2] = bTNode_ref_temp;
	bTNode_AnchorMap_MainHold[regionNum - 2].ID = regionNum - 2;
	int location = 0;
	for (int i = 0; i<regionNum_anchor; i++)
		if (bTNode_AnchorMap[i].hasMerged == false)
		{
			bTNode_AnchorMap_MainHold[location] = bTNode_AnchorMap[i];
			bTNode_AnchorMap_MainHold[location].ID = location;
			location++;
		}
	if (location != regionNum - 2)
	{
		printf("fill location wrong!");
		system("pause");
	}
	//���ref
	BTNode* bTNode_AnchorMap_RefHold = new BTNode[regionNum];
	bTNode_main_temp = bTNode_AnchorMap[cConflictPair.onlyMainNode[0]];
	bTNode_ref_temp = bTNode_AnchorMap[cConflictPair.onlyRefNode[0]];
	for (int i = 1; i<cConflictPair.onlyMainNode.size(); i++)
		bTNode_main_temp = bTNode_main_temp + bTNode_AnchorMap[cConflictPair.onlyMainNode[i]];
	for (int i = 1; i<cConflictPair.onlyRefNode.size(); i++)
		bTNode_ref_temp = bTNode_ref_temp + bTNode_AnchorMap[cConflictPair.onlyRefNode[i]];
	for (int i = 0; i<cConflictPair.iNNode.size(); i++)
		bTNode_ref_temp = bTNode_ref_temp + bTNode_AnchorMap[cConflictPair.iNNode[i]];
	bTNode_AnchorMap_RefHold[regionNum - 1] = bTNode_main_temp;
	bTNode_AnchorMap_RefHold[regionNum - 1].ID = regionNum - 1;
	bTNode_AnchorMap_RefHold[regionNum - 2] = bTNode_ref_temp;
	bTNode_AnchorMap_RefHold[regionNum - 2].ID = regionNum - 2;
	location = 0;
	for (int i = 0; i<regionNum_anchor; i++)
		if (bTNode_AnchorMap[i].hasMerged == false)
		{
			bTNode_AnchorMap_RefHold[location] = bTNode_AnchorMap[i];
			bTNode_AnchorMap_RefHold[location].ID = location;
			location++;
		}
		if (location != regionNum - 2)
		{
			printf("fill location wrong!");
			system("pause");
		}
	//***����Main_hold��Ref_hold��������
	//***1.ֱ�����������ڲ�ͬ����
	//***2.��������ڲ�ͬ�����������������
	//***3.����ȫ��Ӱ��ָ�������������E��
	//***4.ͬ3��Zhang Z��
	//***ע��Main_hold��Ref_hold����ID�Ͷ�Ӧ������ȷ
	//1
	int mode = 4;
	bool flag = true;
	if (mode == 1)
	{
		BTNodeReInitialization(bTNode_AnchorMap_MainHold[regionNum - 2], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_MainHold[regionNum - 1], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_RefHold[regionNum - 2], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_RefHold[regionNum - 1], srimg);
		double ssd_mainhold[2], ssd_refhold[2];
		ssd_mainhold[0] = GetSpectualStandDeviation(bTNode_AnchorMap_MainHold[regionNum - 1], srimg);
		ssd_mainhold[1] = GetSpectualStandDeviation(bTNode_AnchorMap_MainHold[regionNum - 2], srimg);
		ssd_refhold[0] = GetSpectualStandDeviation(bTNode_AnchorMap_RefHold[regionNum - 1], srimg);
		ssd_refhold[1] = GetSpectualStandDeviation(bTNode_AnchorMap_RefHold[regionNum - 2], srimg);
		if (ssd_mainhold[0] + ssd_mainhold[1] <= ssd_refhold[0] + ssd_refhold[1])
			flag = true;
		else
			flag = false;
	}
	if (mode == 3)
	{
		BTNodeReInitialization(bTNode_AnchorMap_MainHold[regionNum - 2], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_MainHold[regionNum - 1], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_RefHold[regionNum - 2], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_RefHold[regionNum - 1], srimg);
		double E_mainHold = Evalution_E(bTNode_AnchorMap_MainHold, regionNum, srimg);
		double E_refHold = Evalution_E(bTNode_AnchorMap_RefHold, regionNum, srimg);
		printf("evalution main ref dif:%lf ,%lf, %lf\n", E_mainHold, E_refHold, E_refHold - E_mainHold);
		if (E_mainHold - E_refHold < T)
			flag = true;
		else
			flag =  false;
	}
	if (mode == 4)
	{
		BTNodeReInitialization(bTNode_AnchorMap_MainHold[regionNum - 2], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_MainHold[regionNum - 1], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_RefHold[regionNum - 2], srimg);
		BTNodeReInitialization(bTNode_AnchorMap_RefHold[regionNum - 1], srimg);
		double Z_mainHold = Evalution_Z(bTNode_AnchorMap_MainHold, regionNum, srimg);
		double Z_refHold = Evalution_Z(bTNode_AnchorMap_RefHold, regionNum, srimg);
		printf("evalution main ref:%lf ,%lf, %lf\n", Z_mainHold, Z_refHold, Z_refHold- Z_mainHold);
		if (Z_mainHold - Z_refHold < T)
			flag = true;
		else
			flag = false;
	}
	//ͳһ�ͷ��ڴ�
	delete[] bTNode_AnchorMap_MainHold;
	delete[] bTNode_AnchorMap_RefHold;
	return flag;
}

int UpdateAnchorMap(vector<BTNode> & bTNode_AnchorMap, int regionNum_anchor,BTNode* bTNode_basic, int regionNum_basic, vector<CConflictPair> & cConflictPair, Mat srimg, vector<BTNode> & new_AnchorMap)
{
	/*
	*����anchor map
	*�����߶ȼ��ѱ�֤��anchor�����򲻻ᷢ����ͻ
	*cConflictPair�еĻ��������ӦbTNode_AnchorMap�е�����
	*��������������ê��
	*/
	bool* mainAnchor = new bool[cConflictPair.size()];
	//��һ����
	for (int i = 0; i<cConflictPair.size(); i++)
	{
		if (MainIsBetterThanRef(bTNode_AnchorMap, regionNum_anchor, cConflictPair[i], srimg) == true)
			mainAnchor[i] = true;
		else
			mainAnchor[i] = false;
	}
	CheckWhichBeAnchored(mainAnchor, cConflictPair.size());
	for (int i = 0; i<regionNum_anchor; i++)
		bTNode_AnchorMap[i].hasMerged = false;
	//����ê��
	int noAnchor = 1; //����Ƿ�������ê��
	for (int i = 0; i<cConflictPair.size(); i++)
	{
		BTNode bTNode_temp;

		////��ê����������
		//for (int j = 0; j<cConflictPair[i].iNNode.size(); j++)
		//{
		//	bTNode_temp = bTNode_temp + bTNode_AnchorMap[cConflictPair[i].iNNode[j]];
		//	bTNode_AnchorMap[cConflictPair[i].iNNode[j]].hasMerged = true;
		//	bTNode_temp.anchored = true;
		//	bTNode_temp.ID = i;
		//	new_AnchorMap.push_back(bTNode_temp);
		//}
		//main������ê��
		if (mainAnchor[i] == true)
		{
			//ֱ����������ê����
			/*bTNode_temp = bTNode_AnchorMap[cConflictPair[i].onlyMainNode[0]];
			bTNode_AnchorMap[cConflictPair[i].onlyMainNode[0]].hasMerged = true;
			for (int j = 1; j<cConflictPair[i].onlyMainNode.size(); j++)
			{
				bTNode_temp = bTNode_temp + bTNode_AnchorMap[cConflictPair[i].onlyMainNode[j]];
				bTNode_AnchorMap[cConflictPair[i].onlyMainNode[j]].hasMerged = true;
			}
			for (int j = 0; j<cConflictPair[i].iNNode.size(); j++)
			{
				bTNode_temp = bTNode_temp + bTNode_AnchorMap[cConflictPair[i].iNNode[j]];
				bTNode_AnchorMap[cConflictPair[i].iNNode[j]].hasMerged = true;
			}
			bTNode_temp.anchored = true;
			bTNode_temp.ID = i;
			new_AnchorMap.push_back(bTNode_temp);*/
		}
		//main����ʴ��ê��
		else
		{
			noAnchor = 0;
			bTNode_temp = bTNode_AnchorMap[cConflictPair[i].onlyMainNode[0]];
			bTNode_AnchorMap[cConflictPair[i].onlyMainNode[0]].hasMerged = true;
			for (int j = 1; j<cConflictPair[i].onlyMainNode.size(); j++)
			{
				bTNode_temp = bTNode_temp + bTNode_AnchorMap[cConflictPair[i].onlyMainNode[j]];
				bTNode_AnchorMap[cConflictPair[i].onlyMainNode[j]].hasMerged = true;
			}
			bTNode_temp.anchored = true;
			bTNode_temp.ID = i;
			new_AnchorMap.push_back(bTNode_temp);
		}
	}
	//δ����ϲ������
	for (int i = 0; i<regionNum_anchor; i++)
		if (bTNode_AnchorMap[i].hasMerged == false)
			new_AnchorMap.push_back(bTNode_AnchorMap[i]);
	return noAnchor;
}

int MergeErrorDetention_Main(BTNode* bTNode_basic, int regionNum_basic, BTNode* bTNode_main, BTNode* bTNode_ref, int regionNum_main, int regionNum_ref, Mat srimg, vector<BTNode> & bTNode_newAnchorMap)
{
	/*
	*�ϲ������� �ӿں���
	*/
	//��ʼ��Anchor Map
	CheckMergeErrorDetentionStart();
	vector<BTNode> bTNode_AnchorMap;
	for (int i = 0; i<regionNum_basic; i++)
	{
		//bTNode_basic[i].anchored = false;      
		bTNode_basic[i].hasMerged = false;
		bTNode_AnchorMap.push_back(bTNode_basic[i]);
	}
	vector<CConflictPair> cConflictPair;
	int regionNum_AnchorMap = bTNode_AnchorMap.size();
	//ƥ�䲢����ͻ����
	FindRegionConflict(bTNode_AnchorMap, regionNum_AnchorMap, bTNode_main, bTNode_ref, regionNum_main, regionNum_ref, cConflictPair, srimg);
	CheckConflictRegionNum(cConflictPair.size());
	if (cConflictPair.size() == 0)
	{
		printf("no conflict!\n");
		return 0;
	}
	int noAnchor = 1;
	noAnchor = UpdateAnchorMap(bTNode_AnchorMap, regionNum_AnchorMap, bTNode_basic, regionNum_basic, cConflictPair, srimg, bTNode_newAnchorMap);
	if (noAnchor == 1)
	{
		printf("no segment anchor!\n");
		return 0;
	}
	return 1;
}
#endif