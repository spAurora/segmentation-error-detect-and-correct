//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////
#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "ClassAndCheck.h"

void RecursivelySetNowNodeValue(int* labels, int & setValue, BTNode* bTNode)
{
	/*
	*�ݹ���������ֵ
	*/
	if (bTNode->left == NULL && bTNode->right == NULL)
	{
		vector<int>::iterator it;
		for (it = bTNode->pixelLocation.begin(); it != bTNode->pixelLocation.end(); it++)
			labels[*it] = setValue;
	}
	if (bTNode->left != NULL)
		RecursivelySetNowNodeValue(labels, setValue, bTNode->left);
	if (bTNode->right != NULL)
		RecursivelySetNowNodeValue(labels, setValue, bTNode->right);
}

void RecursivelySetAllNodeValue(int*labels, BTNode* bTNode, int regionNum, int & setValue)
{
	/*
	*�ݹ�����������
	*/
	if (bTNode->left == NULL && bTNode->right == NULL)
	{
		setValue++;
		RecursivelySetNowNodeValue(labels, setValue, bTNode);
	}
	else if(bTNode->left->highLevelNode == false && bTNode->right->highLevelNode == false)
	{
		setValue++;
		RecursivelySetNowNodeValue(labels, setValue, bTNode);
	}
	if (bTNode->left != NULL)
		if (bTNode->left->highLevelNode == true)
			RecursivelySetAllNodeValue(labels, bTNode->left, regionNum, setValue);
	if (bTNode->right != NULL)
		if (bTNode->right->highLevelNode == true)
			RecursivelySetAllNodeValue(labels, bTNode->right, regionNum, setValue);
}

void SetAllNodeValueByAbsoluteLocation(int* labels, BTNode* bTnode, int regionNum, int absoluteLocation)
{
	/*
	*���ھ��Ժϲ�λ�ü���ϲ����
	*����ͨͼ����ʽ����
	*/
	//ˢ��ԭ����
	for (int i = 0; i<2*regionNum-1; i++)
		bTnode[i].highLevelNode = true;
	if (absoluteLocation<0 || absoluteLocation>regionNum-1)
	{
		printf("Location Wrong!");
		system("pause");
		return;
	}
	//����ǰ�κϲ�����Ϊ���ɼ�
	for (int i = regionNum; i<regionNum+absoluteLocation; i++)
	{
		bTnode[bTnode[i].left->ID].highLevelNode = false;
		bTnode[bTnode[i].right->ID].highLevelNode = false;
	}
	int setValue = -1;
	RecursivelySetAllNodeValue(labels, &bTnode[2*regionNum-2], regionNum, setValue);
}

void SetAllNodeValueByRelativeLocation(int* labels, BTNode* bTnode, int regionNum, double relativeLocation)
{
	/*
	*������Ժϲ�λ�ü���ϲ����
	*����ͨͼ����ʽ����
	*/
	//ˢ��ԭ����
	for (int i = 0; i<2*regionNum-1; i++)
		bTnode[i].highLevelNode = true;
	if (relativeLocation<0 || relativeLocation>1)
	{
		printf("Location Wrong!");
		system("pause");
		return;
	}
	//���λ��ת��Ϊ����λ��
	int absoluteLocation = (int)(regionNum-1)*relativeLocation;
	//����ǰ�κϲ�����Ϊ���ɼ�
	for (int i = regionNum; i<regionNum+absoluteLocation; i++)
	{
		bTnode[bTnode[i].left->ID].highLevelNode = false;
		bTnode[bTnode[i].right->ID].highLevelNode = false;
	}
	int setValue = -1;
	RecursivelySetAllNodeValue(labels, &bTnode[2*regionNum-2], regionNum, setValue);
}

void OutputSegmentResult(int* labels, Mat & srimg, char output_Path[], char output_ImgName[])
{
	/*
	*����labels���Ӱ��ָ���
	*/
	int width = srimg.cols;
	int height = srimg.rows;
	Mat	imgOut = srimg.clone();
	for (int i = 1; i<height-1; i++)
		for (int j = 1;j<width-1; j++)
		{
			if (labels[i*width + j] != labels[(i+1)*width +j] || labels[i*width + j] != labels[i*width +j+1])
			{
				imgOut.data[(i*width + j)*3] = 0;
				imgOut.data[(i*width + j)*3 + 1] = 0;
				imgOut.data[(i*width + j)*3 + 2] = 255;
			}
		}
	char output_Path_new[1024], output_Path_txt[1024];
	strcpy(output_Path_new, output_Path);
	strcat(output_Path_new, output_ImgName);
	imwrite(output_Path_new, imgOut);

	strcpy(output_Path_txt, output_Path);
	strcat(output_Path_txt, output_ImgName);
	strcat(output_Path_txt, ".txt");
	FILE *fp;
	if((fp = fopen(output_Path_txt, "w+")) == NULL)
	{
		printf("��outlabels�ļ�ʧ��\n");
		system("pause");
		exit(-1);
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
			fprintf(fp, "%d\t", labels[i*width + j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void OutputAnchorMap(int* labels, Mat & srimg, char output_Path[], char output_ImgName[], vector<BTNode> & bTNode_AnchorMap)
{
	/*
	*ר��������ʾanchormap��ê�������Ե��һ����ɫ
	*/
	int width = srimg.cols;
	int height = srimg.rows;
	Mat	imgOut = srimg.clone();
	for (int i = 1; i<height-1; i++)
		for (int j = 1;j<width-1; j++)
		{
			if (labels[i*width + j] != labels[(i+1)*width +j] || labels[i*width + j] != labels[i*width +j+ 1])
			{
				if (bTNode_AnchorMap[labels[i*width + j]].anchored == true || bTNode_AnchorMap[labels[(i+1)*width +j]].anchored == true || bTNode_AnchorMap[labels[(i+1)*width +j +1]].anchored == true)
				{
					//ê�������Ե��ʾ��ɫ
					imgOut.data[(i*width + j)*3] = 232;
					imgOut.data[(i*width + j)*3 + 1] = 162;
					imgOut.data[(i*width + j)*3 + 2] = 0;
				}
				else
				{
					//���������Ե��ʾ��ɫ
					imgOut.data[(i*width + j)*3] = 0;
					imgOut.data[(i*width + j)*3 + 1] = 0;
					imgOut.data[(i*width + j)*3 + 2] = 255;
				}
			}
		}
	char output_Path_new[1024];
	strcpy(output_Path_new, output_Path);
	strcat(output_Path_new, output_ImgName);
	imwrite(output_Path_new, imgOut);
}

void OutputAnchorMap_fill(int* labels, Mat & srimg, char output_Path[], char output_ImgName[], vector<BTNode> & bTNode_AnchorMap)
{
	/*
	*ר��������ʾanchormap��ê���������Ϊ��һ����ɫ
	*/
	int width = srimg.cols;
	int height = srimg.rows;
	Mat	imgOut = srimg.clone();
	for (int i = 1; i<height-1; i++)
		for (int j = 1;j<width-1; j++)
		{
			if (labels[i*width + j] != labels[(i+1)*width +j] || labels[i*width + j] != labels[i*width +j+ 1])
			{
					//���������Ե��ʾ��ɫ
					imgOut.data[(i*width + j)*3] = 0;
					imgOut.data[(i*width + j)*3 + 1] = 0;
					imgOut.data[(i*width + j)*3 + 2] = 255;
			}
		}
	//ê���������
	for (int i = 0; i<height; i++)
		for (int j = 0; j<width; j++)
		{
			if (bTNode_AnchorMap[labels[i*width + j]].anchored == true)
			{
				imgOut.data[(i*width + j)*3] = 232;
				imgOut.data[(i*width + j)*3 + 1] = 162;
				imgOut.data[(i*width + j)*3 + 2] = 0;
			}
		}
	char output_Path_new[1024];
	strcpy(output_Path_new, output_Path);
	strcat(output_Path_new, output_ImgName);
	imwrite(output_Path_new, imgOut);
}

void SetLabelsValueByBTNodeList(vector<BTNode> & bTNode, int* labels)
{
	/*
	*����bTNode List����labels
	*/
	for (int i = 0; i<bTNode.size(); i++)
		for (int j = 0; j<bTNode[i].pixelLocation.size(); j++)
			labels[bTNode[i].pixelLocation[j]] = i; 
}

#endif