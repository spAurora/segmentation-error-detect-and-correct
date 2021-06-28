//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////
#ifndef _EVALUTION_H_
#define _EVALUTION_H_

#include "ClassAndCheck.h"
#include "ScaleSets.h"
#include "Output.h"

using namespace std;
using namespace cv;

#define NORMALIZATIONTERM 10

double GetBrightness(BTNode bTnode)
{
	/*
	*��������ƽ������
	*/
	return (bTnode.avgB + bTnode.avgG + bTnode.avgR)/3;
}

double GetSpectualStandDeviation(BTNode bTNode, Mat srimg)
{
	/*
	*���ع��ױ�׼��
	*/
	double tempBsqr = 0, tempGsqr = 0, tempRsqr = 0;
	double tempSqr = 0, sumSqr = 0;
	for (int j = 0; j<bTNode.pixelLocation.size(); j++)
	{
		tempBsqr = (srimg.data[bTNode.pixelLocation[j]*3] - bTNode.avgB)*(srimg.data[bTNode.pixelLocation[j]*3] - bTNode.avgB);
		tempGsqr = (srimg.data[bTNode.pixelLocation[j]*3+1] - bTNode.avgG)*(srimg.data[bTNode.pixelLocation[j]*3+1] - bTNode.avgG);
		tempRsqr = (srimg.data[bTNode.pixelLocation[j]*3+2] - bTNode.avgR)*(srimg.data[bTNode.pixelLocation[j]*3+2] - bTNode.avgR);
		tempSqr = (tempBsqr + tempGsqr + tempRsqr)/3;
		sumSqr += tempSqr;
	}
	return sqrt(sumSqr/bTNode.area);
}

void CalculateVKMI(int* labels, BTNode* bTNode, int regionNum, Mat srimg, double & VK, double & MI, int location)
{
	/*
	*ͬ������VK��MI
	*/
	int regionNum_new = 0;
	int height = srimg.rows;
	int width = srimg.cols;
	for (int i = 0; i<height; i++)
		for (int j = 0; j<width; j++)
			if (labels[i*width + j] > regionNum_new)
				regionNum_new = labels[i*width + j];
	regionNum_new++;
	for (int i = 0; i<2*regionNum-1; i++)
		bTNode[i].highLevelNode = true;
	for (int i = regionNum; i < regionNum + location; i++)
	{
		bTNode[bTNode[i].left->ID].highLevelNode = false;
		bTNode[bTNode[i].right->ID].highLevelNode = false;
	}
	//����ǰ����������֯
	//ע�����п�����ָ�룬�����ٵ���
	BTNode* bTNode_new = new BTNode[regionNum_new];
	int fill_location = 0;
	for (int i = 0; i<regionNum + location; i++)
		if (bTNode[i].highLevelNode == true)
		{
			bTNode_new[fill_location] = bTNode[i];
			fill_location++;
		}
	if (fill_location != regionNum_new)
	{
		printf("fill location wrong!\n");
		system("pause");
	}
	//����VK
	double sumaq = 0, tempaq = 0, spectualStandDeviation = 0;
	for (int i = 0; i<regionNum + location; i++)
	{
		double tempBsqr = 0, tempGsqr = 0, tempRsqr = 0;
		double tempSqr = 0, sumSqr = 0;
		if (bTNode[i].highLevelNode == true)
		{
			for (int j = 0; j<bTNode[i].pixelLocation.size(); j++)
			{
				tempBsqr = (srimg.data[bTNode[i].pixelLocation[j]*3] - bTNode[i].avgB)*(srimg.data[bTNode[i].pixelLocation[j]*3] - bTNode[i].avgB);
				tempGsqr = (srimg.data[bTNode[i].pixelLocation[j]*3+1] - bTNode[i].avgG)*(srimg.data[bTNode[i].pixelLocation[j]*3+1] - bTNode[i].avgG);
				tempRsqr = (srimg.data[bTNode[i].pixelLocation[j]*3+2] - bTNode[i].avgR)*(srimg.data[bTNode[i].pixelLocation[j]*3+2] - bTNode[i].avgR);
				tempSqr = (tempBsqr + tempGsqr + tempRsqr)/3;
				sumSqr += tempSqr;
			}
			spectualStandDeviation = sqrt(sumSqr/bTNode[i].area);
			//������Թ��ױ�׼��
			tempaq = bTNode[i].area*(spectualStandDeviation); 
			sumaq += tempaq;
		}
	}
	VK = sumaq / (width * height);
	//����MI
	double x_meancolor = 0;
	double sumB = 0, sumG = 0, sumR = 0;
	for (int i = 0;i<height;i++)
		for (int j = 0; j<width; j++)
		{
			sumB += srimg.data[(i*width+j)*3];
			sumG += srimg.data[(i*width+j)*3+1];
			sumR += srimg.data[(i*width+j)*3+2];
		}
		x_meancolor = (sumB + sumG +sumR)/(3*(height*width));
	//�½�����ͼ
	ArrayHeadGraphNode* head = new ArrayHeadGraphNode[regionNum_new];
	CreateToplogicalGraph(labels, head, regionNum_new, width, height);
	double sumXijmean = 0, sumXii = 0;
	int sumWij = 0;
	for (int i = 0; i<regionNum_new; i++)
	{
		vector<GraphNode>::iterator it;
		for (it = head[i].pGraphNodeList.begin(); it != head[i].pGraphNodeList.end(); it++)
		{
			//�ڽӼ���
			sumWij++; 
			sumXijmean += (GetBrightness(bTNode_new[i]) - x_meancolor) * (GetBrightness(bTNode_new[it->ID]) - x_meancolor);
		}
	}
	//�ظ�ͳ����2��
	sumXijmean /= 2;
	sumWij /= 2;
	for (int i = 0; i<regionNum_new; i++)
		sumXii += (GetBrightness(bTNode_new[i]) - x_meancolor)*(GetBrightness(bTNode_new[i]) - x_meancolor);
	MI = (regionNum_new*sumXijmean) / (sumWij*sumXii);
}

double CalculateVK(BTNode* bTNode, int regionNum, Mat srimg)
{
	/*
	*����VK
	*/
	double sumaq = 0, tempaq = 0, spectualStandDeviation = 0;
	int width =srimg.cols;
	int height = srimg.rows;
	for (int i = 0; i<regionNum; i++)
	{
		double tempBsqr = 0, tempGsqr = 0, tempRsqr = 0;
		double tempSqr = 0, sumSqr = 0;
		for (int j = 0; j<bTNode[i].pixelLocation.size(); j++)
		{
			tempBsqr = (srimg.data[bTNode[i].pixelLocation[j]*3] - bTNode[i].avgB)*(srimg.data[bTNode[i].pixelLocation[j]*3] - bTNode[i].avgB);
			tempGsqr = (srimg.data[bTNode[i].pixelLocation[j]*3+1] - bTNode[i].avgG)*(srimg.data[bTNode[i].pixelLocation[j]*3+1] - bTNode[i].avgG);
			tempRsqr = (srimg.data[bTNode[i].pixelLocation[j]*3+2] - bTNode[i].avgR)*(srimg.data[bTNode[i].pixelLocation[j]*3+2] - bTNode[i].avgR);
			tempSqr = (tempBsqr + tempGsqr + tempRsqr)/3;
			sumSqr += tempSqr;
		}
		spectualStandDeviation = sqrt(sumSqr/bTNode[i].area);
			//������Թ��ױ�׼��
		tempaq = bTNode[i].area*(spectualStandDeviation); 
		sumaq += tempaq;
	}
	return sumaq / (width * height);
}

double CalculateMI(BTNode* bTNode, int regionNum, Mat srimg)
{
	/*
	*����MI
	*/
	double x_meancolor = 0;
	double sumB = 0, sumG = 0, sumR = 0;
	int width =srimg.cols;
	int height = srimg.rows;
	for (int i = 0;i<height;i++)
		for (int j = 0; j<width; j++)
		{
			sumB += srimg.data[(i*width+j)*3];
			sumG += srimg.data[(i*width+j)*3+1];
			sumR += srimg.data[(i*width+j)*3+2];
		}
		x_meancolor = (sumB + sumG +sumR)/(3*(height*width));
	//�½�labelsͼ��
	int* labels = new int[width*height];
	for (int i = 0; i<regionNum; i++)
		for (int j = 0; j<bTNode[i].pixelLocation.size(); j++)
			labels[bTNode[i].pixelLocation[j]] = i;
	//�½�����ͼ
	ArrayHeadGraphNode* head = new ArrayHeadGraphNode[regionNum];
	CreateToplogicalGraph(labels, head, regionNum, width, height);
	double sumXijmean = 0, sumXii = 0;
	int sumWij = 0;
	for (int i = 0; i<regionNum; i++)
	{
		vector<GraphNode>::iterator it;
		for (it = head[i].pGraphNodeList.begin(); it != head[i].pGraphNodeList.end(); it++)
		{
			//�ڽӼ���
			sumWij++; 
			sumXijmean += (GetBrightness(bTNode[i]) - x_meancolor) * (GetBrightness(bTNode[it->ID]) - x_meancolor);
		}
	}
	//�ظ�ͳ����2��
	sumXijmean /= 2;
	sumWij /= 2;
	for (int i = 0; i<regionNum; i++)
		sumXii += (GetBrightness(bTNode[i]) - x_meancolor)*(GetBrightness(bTNode[i]) - x_meancolor);
	return (regionNum*sumXijmean) / (sumWij*sumXii);
}

int CalculateBestScaleLocation(BTNode* bTNode, int regionNum, Mat srimg, double Q, BTNode* bTNode_LevelOne, int regionNum_LevelOne)
{
	/*
	*���ּ������ų߶�λ��
	*�Ժϲ����������������Ϊ���պϲ�
	*/
	//����VK MI first
	double VK_first, MI_first, VK_final, MI_final;
	VK_first = CalculateVK(bTNode_LevelOne, regionNum_LevelOne, srimg);
	MI_first = CalculateMI(bTNode_LevelOne, regionNum_LevelOne, srimg);
	//����bTNode�ϲ�λ��
	for (int i = 0; i<2*regionNum-1; i++)
		bTNode[i].highLevelNode = true;
	//�Ե����ڶ���λ��Ϊ���պϲ�λ��
	for (int i = regionNum; i < 2*regionNum - 2; i++)
	{
		bTNode[bTNode[i].left->ID].highLevelNode = false;
		bTNode[bTNode[i].right->ID].highLevelNode = false;
	}
	//��ʱȫ��ʣ��������
	BTNode* bTNode_final = new BTNode[2];
	int fill_location = 0;
	for (int i = 0; i < 2*regionNum - 2; i++)
		if (bTNode[i].highLevelNode == true)
		{
			bTNode_final[fill_location] = bTNode[i];
			fill_location++;
		}
	if (fill_location != 2)
	{
		printf("fill location wrong!\n");
		system("pause");
	}
	
	//����VK MI final
	VK_final = CalculateVK(bTNode_final, 2, srimg);
	MI_final = CalculateMI(bTNode_final, 2, srimg);
	CheckVKMIFirstFinal(VK_first, VK_final, MI_first, MI_final);

	int low = 0;
	int high = regionNum - 2;
	int mid = (low + high)/2;
	int width = srimg.cols;
	int height = srimg.rows;
	double VK, MI, PU, PO;
	while (low < high)
	{
		mid = (low + high)/2;
		//����bTNode�ϲ�λ��
		for (int i = 0; i<2*regionNum-1; i++)
			bTNode[i].highLevelNode = true;
		//��midΪ���պϲ�λ��
		for (int i = regionNum; i < regionNum + mid; i++)
		{
			bTNode[bTNode[i].left->ID].highLevelNode = false;
			bTNode[bTNode[i].right->ID].highLevelNode = false;
		}
		//����һ�κϲ�����������-1
		BTNode* bTNode_now = new BTNode[regionNum - mid];
		int fill_location_main = 0;
		for (int i = 0; i < regionNum + mid; i++)
			if (bTNode[i].highLevelNode == true)
			{
				bTNode_now[fill_location_main] = bTNode[i];
				fill_location_main++;
			}
		if (fill_location_main != regionNum - mid)
		{
			printf("fill location main wrong!\n");
			system("pause");
		}
		//����VK MI
		VK = CalculateVK(bTNode_now, regionNum - mid, srimg);
		MI = CalculateMI(bTNode_now, regionNum - mid, srimg);
		//��һ��
		PU = (VK - VK_first) / (VK_final - VK_first);
		PO = (MI - MI_final) / (MI_first - MI_final);
		if (PU*Q == PO)
			break;
		if (PU*Q < PO)
			low = mid + 1;
		else if (PU*Q > PO)
			high = mid - 1;
		printf("Q*PU PO: %lf, %lf\n", Q*PU, PO);
		//���뼰ʱ�ͷ�
		delete[] bTNode_now;
	} 
	printf("low mid high:%d %d %d\n", low, mid, high);
	return low;
}

double GetColorVectorError_e(BTNode bTNode, Mat srimg)
{
	/*
	*��������Zhang�����������е�e
	*/
	double tempBsqr = 0, tempGsqr = 0, tempRsqr = 0;
	double tempSqr = 0, sumSqr = 0;
	for (int j = 0; j<bTNode.pixelLocation.size(); j++)
	{
		tempBsqr = (srimg.data[bTNode.pixelLocation[j]*3] - bTNode.avgB)*(srimg.data[bTNode.pixelLocation[j]*3] - bTNode.avgB);
		tempGsqr = (srimg.data[bTNode.pixelLocation[j]*3+1] - bTNode.avgG)*(srimg.data[bTNode.pixelLocation[j]*3+1] - bTNode.avgG);
		tempRsqr = (srimg.data[bTNode.pixelLocation[j]*3+2] - bTNode.avgR)*(srimg.data[bTNode.pixelLocation[j]*3+2] - bTNode.avgR);
		sumSqr += sqrt(tempBsqr + tempGsqr + tempRsqr);
	}
	sumSqr /= bTNode.area;
	return sumSqr;
}
double CalculateLambada(BTNode* bTNode, int regionNum, Mat srimg)
{	
	/*
	*��ѧ�� �޼ල�ָ�����Z
	*ϵ��lambada����
	*�����ʼ�ָ���
	*/
	double e;
	BTNode bTNode_temp;
	int width = srimg.cols;
	int height = srimg.rows;
	int S = width*height;
	//bTNode_temp��Ϊ�ϲ���������
	bTNode_temp.area = width*height;
	for (int i = 0; i<regionNum; i++)
		for (int j = 0; j<bTNode[i].pixelLocation.size(); j++)
			bTNode_temp.pixelLocation.push_back(bTNode[i].pixelLocation[j]);
	double tempB = 0, tempG = 0, tempR = 0;
	for (int i = 0; i<bTNode_temp.pixelLocation.size(); i++)
	{
		tempB += srimg.data[bTNode_temp.pixelLocation[i] * 3];
		tempG += srimg.data[bTNode_temp.pixelLocation[i] * 3 + 1];
		tempR += srimg.data[bTNode_temp.pixelLocation[i] * 3 + 2];
	}
	bTNode_temp.avgB = tempB/S;
	bTNode_temp.avgG = tempG/S;
	bTNode_temp.avgR = tempR/S;
	e = GetColorVectorError_e(bTNode_temp, srimg);
	//printf("bTNode_temp BGR: %lf, %lf, %lf, ... e:%lf\n",bTNode_temp.avgB, bTNode_temp.avgG, bTNode_temp.avgR, e);
	double Tmin = (1/(double)(NORMALIZATIONTERM)) * e / (1 + log((double)bTNode_temp.area));
	
	int R = regionNum;
	double totol_term3 = 0;
	//��������ۼ��� log��eΪ��
	for (int i = 0; i<regionNum; i++)
	{
		e = GetColorVectorError_e(bTNode[i], srimg);
		totol_term3 += e / (1 + log((double)bTNode[i].area));
	}
	double Tmax;
	Tmax = (1/(double)(NORMALIZATIONTERM)) * sqrt((double)regionNum) * totol_term3;
	
	double Dmin = 0;
	double Dmax;
	double avgB = 0, avgG = 0, avgR = 0;
	for (int i = 0; i<height; i++)
		for (int j = 0; j<width; j++)
		{
			avgB += srimg.data[3*(i*width+height)];
			avgG += srimg.data[3*(i*width+height) + 1];
			avgR += srimg.data[3*(i*width+height) + 2];
		}
	avgB /= (width*height);
	avgG /= (width*height);
	avgR /= (width*height);
	double D_term1 = 0;
	for (int i =0; i<regionNum; i++)
	{
		D_term1 += (bTNode[i].avgB - avgB) * (bTNode[i].avgB - avgB);
		D_term1 += (bTNode[i].avgG - avgG) * (bTNode[i].avgG - avgG);
		D_term1 += (bTNode[i].avgR - avgR) * (bTNode[i].avgR - avgR);
	}
	D_term1 /= regionNum;
	Dmax = D_term1 / sqrt((double)regionNum);
	printf("Tmin:%lf, Tmax:%lf, Dmin:%lf, Dmax:%lf\n", Tmin, Tmax, Dmin, Dmax);
	return (Tmax - Tmin) / (Dmax - Dmin);
}

double Evalution_Z(BTNode* bTNode, int regionNum, Mat srimg)
{
	/*
	*��ѧ�� �޼ල�ָ�����Z ����
	*/
	int width = srimg.cols;
	int height = srimg.rows;
	int S = width*height;
	int R = regionNum;
	double totol_term3 = 0;
	//��������ۼ��� log��eΪ��
	for (int i = 0; i<regionNum; i++)
	{
		double e = 0;
		e = GetColorVectorError_e(bTNode[i], srimg);
		totol_term3 += e / (1 + log((double)bTNode[i].area));
	}
	double T;
	//����T
	T = (1/(double)(NORMALIZATIONTERM)) * sqrt((double)regionNum) * totol_term3;
	
	double avgB = 0, avgG = 0, avgR = 0;
	for (int i = 0; i<height; i++)
		for (int j = 0; j<width; j++)
		{
			avgB += srimg.data[3*(i*width+height)];
			avgG += srimg.data[3*(i*width+height) + 1];
			avgR += srimg.data[3*(i*width+height) + 2];
		}
	avgB /= (width*height);
	avgG /= (width*height);
	avgR /= (width*height);
	double D_term1 = 0;
	for (int i =0; i<regionNum; i++)
	{
		D_term1 += (bTNode[i].avgB - avgB) * (bTNode[i].avgB - avgB);
		D_term1 += (bTNode[i].avgG - avgG) * (bTNode[i].avgG - avgG);
		D_term1 += (bTNode[i].avgR - avgR) * (bTNode[i].avgR - avgR);
	}
	//����
	D_term1 /= regionNum;    
	double D;
	//����D
	D = D_term1 / sqrt((double)regionNum);
	//����Z
	double Z;
	double lambada;
	lambada = CalculateLambada(bTNode, regionNum, srimg);
	Z = T + lambada*D;
	printf("T:%lf, D:%lf, l*D:%lf\n", T, D, lambada*D);
	return Z;
}

void CalculateL(int* L, BTNode bTNode, int band, Mat srimg)
{
	/*
	*ͳ��L�����ڼ���Hv
	*/
	for (int i = 0; i<bTNode.pixelLocation.size(); i++)
		L[srimg.data[bTNode.pixelLocation[i]*3 + band]]++;
}

double Evalution_E(BTNode* bTNode, int regionNum, Mat srimg)
{	
	/*
	*Hui Zhang
	*�����ص�Ӱ��ָ�����
	*/
	int SI = srimg.rows * srimg.cols;
	double Hr_B = 0, Hr_G = 0, Hr_R = 0;
	for (int i =0; i<regionNum; i++)
	{
		int* L_B = new int[256] ();
		int* L_G = new int[256] ();
		int* L_R = new int[256] ();
		CalculateL(L_B, bTNode[i], 0, srimg);
		CalculateL(L_G, bTNode[i], 1, srimg);
		CalculateL(L_R, bTNode[i], 2, srimg);
		double Hv_B = 0, Hv_G = 0, Hv_R = 0;
		for (int j = 0; j<256; j++)
		{
			if (L_B[j] != 0)
				Hv_B += -(double)L_B[j]/bTNode[i].area * log((double)L_B[j]/bTNode[i].area);
			if (L_G[j] != 0)
				Hv_G += -(double)L_G[j]/bTNode[i].area * log((double)L_G[j]/bTNode[i].area);
			if (L_R[j] != 0)
				Hv_R += -(double)L_R[j]/bTNode[i].area * log((double)L_R[j]/bTNode[i].area);
		}
		Hr_B += (double)bTNode[i].area/SI * Hv_B;
		Hr_G += (double)bTNode[i].area/SI * Hv_G;
		Hr_R += (double)bTNode[i].area/SI * Hv_R;
		delete[] L_B;
		delete[] L_G;
		delete[] L_R;
	}
	double Hv;
	Hv = (Hr_B + Hr_G + Hr_R) / 3;


	//����Hl
	double Hl = 0;
	for (int i = 0; i<regionNum; i++)
		Hl += -(double)bTNode[i].area/SI * log((double)bTNode[i].area/SI);

	//printf("Hv Hl:%lf, %lf\n", Hv, Hl);
	//����E
	return Hl + Hv;
}

double Evalution_E_Main(BTNode* bTNode, int regionNum, Mat srimg, int location)
{
	/*
	*Hui Zhang
	*�����ص�Ӱ��ָ�����
	*�ӿں���
	*/
	for (int i = 0; i<2*regionNum-1; i++)
		bTNode[i].highLevelNode = true;
	for (int i = regionNum; i < regionNum + location; i++)
	{
		bTNode[bTNode[i].left->ID].highLevelNode = false;
		bTNode[bTNode[i].right->ID].highLevelNode = false;
	}
	BTNode* bTNode_now = new BTNode[regionNum - location];
	int fill_location_main = 0;
	for (int i = 0; i < regionNum + location; i++)
		if (bTNode[i].highLevelNode == true)
		{
			bTNode_now[fill_location_main] = bTNode[i];
			fill_location_main++;
		}
	if (fill_location_main != regionNum - location)
		{
			printf("fill location main wrong!\n");
			system("pause");
		}
	double E;
	E = Evalution_E(bTNode_now, regionNum - location, srimg);
	return E;
}

#endif