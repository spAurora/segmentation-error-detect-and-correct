//////////////////////////////////////////////////////////////////////////
//	code by wHy
//  Aerospace Information Research Institute, Chinese Academy of Sciences
//	751984964@qq.com
//////////////////////////////////////////////////////////////////////////
#ifndef _CRITERION_H_
#define _CRITERION_H_

#include "ClassAndCheck.h"

#define min(a, b) ((a)<(b)?(a):(b))
#define max(a,b,c) (((((a)>(b)) ? (a):(b))>(c)) ? (((a)>(b)) ? (a):(b)):(c))
#define MINC 9999999999999

using namespace std;
using namespace cv;

double CalcuteCriterion_AverageSpectralDifference(BTNode* t1, BTNode* t2, Mat srimg)
{
	/*
	*������׼��_1
	*ƽ�����ײ�
	*/
	return sqrt( ((t1->avgB-t2->avgB)*(t1->avgB-t2->avgB) + (t1->avgG-t2->avgG)*(t1->avgG-t2->avgG) + (t1->avgR-t2->avgR)*(t1->avgR-t2->avgR)) / 3);
}

double CalcuteCriterion_MaxSpectralDifference(BTNode* t1, BTNode* t2, Mat srimg)
{
	/*
	*������׼��_2
	*�����ײ�
	*/
	return max(sqrt((t1->avgB-t2->avgB)*(t1->avgB-t2->avgB)), sqrt((t1->avgG-t2->avgG)*(t1->avgG-t2->avgG)), sqrt((t1->avgR-t2->avgR)*(t1->avgR-t2->avgR)));
}

void CalcutePerimeter(int location, int* label, bool* visited, int & l, int width, int height)
{
	/*
	*�ݹ�ͳ���ܳ�
	*4��ͨģʽ
	*/
	visited[location] = true;
	//�ж��Ƿ�Ӧͳ��Ϊ�ܳ�����
	int l_height = location / width;
	int l_width = location % width;
	//��Եֱ��ͳ��Ϊ�ܳ�����
	if (l_height == 0 || l_height == height-1 || l_width == 0 || l_width == width-1)
	{
		l++;
	}
	else
	{
		//������������ĸ�λ�����ص����
		if (label[l_height * width + l_width - 1] == 1 && label[l_height * width + l_width + 1] && label[(l_height-1) * width + l_width] && label[(l_height+1) * width + l_width])
		{
			//�ڲ���������
		}
		else
			l++;
	}
	//���ĸ�����ݹ����
	if(l_height != 0)
		if(visited[(l_height-1) * width + l_width] == false)
			CalcutePerimeter((l_height-1) * width + l_width, label, visited, l, width, height);
	if(l_height != height-1)
		if(visited[(l_height+1) * width + l_width] == false)
			CalcutePerimeter((l_height+1) * width + l_width, label, visited, l, width, height);
	if(l_width != 0)
		if(visited[l_height*width + l_width - 1] == false)
			CalcutePerimeter(l_height * width + l_width - 1, label, visited, l, width, height);
	if(l_width != width-1)
		if(visited[l_height*width + l_width + 1] == false)
			CalcutePerimeter(l_height*width + l_width + 1, label, visited, l, width, height);
}

void CalculateSpectralStandardDeviation(BTNode* t, double (&SD)[3], Mat srimg)
{
	/*
	*������ױ�׼��
	*�����¼3�����ι��ױ�׼������������
	*/
	double tempBsqr = 0, tempGsqr = 0, tempRsqr = 0;
	double sumBqr = 0, sumGqr = 0, sumRqr = 0;
	for (int i = 0; i<t->pixelLocation.size(); i++)
	{
		tempBsqr = (srimg.data[t->pixelLocation[i]*3] - t->avgB)*(srimg.data[t->pixelLocation[i]*3] - t->avgB);
		tempGsqr = (srimg.data[t->pixelLocation[i]*3+1] - t->avgG)*(srimg.data[t->pixelLocation[i]*3+1] - t->avgG);
		tempRsqr = (srimg.data[t->pixelLocation[i]*3+2] - t->avgR)*(srimg.data[t->pixelLocation[i]*3+2] - t->avgR);

		sumBqr += tempBsqr;
		sumGqr += tempGsqr;
		sumRqr += tempRsqr;
	}
	SD[0] = sqrt(sumBqr / t->area);
	SD[1] = sqrt(sumGqr / t->area);
	SD[2] = sqrt(sumRqr / t->area);
}

void CalculateSpectralStandardDeviation_Fast(BTNode* t, double (&SD)[3], Mat srimg)
{

}

int CalculateBoundingBoxPerimete(BTNode* t, Mat srimg)
{
	/*
	*�����Χ�е��ܳ�
	*/
	int minW = MINC, minH = MINC, maxW = 0, maxH = 0;
	int width = srimg.cols;
	int height = srimg.rows;
	int tempW, tempH;
	for (int i = 0; i<t->pixelLocation.size(); i++)
	{
		tempH = t->pixelLocation[i] / width;
		tempW = t->pixelLocation[i] % width;
		if (tempH < minH)
			minH = tempH;
		if (tempW < minW)
			minW = tempW;
		if (tempH > maxH)
			maxH = tempH;
		if (tempW > maxW)
			maxW = tempW;
	}
	return 2*(maxW - minW + maxH - minH);
}

int MyfindContours(const Mat& src, vector<vector<Point>>& contours, vector<Vec4i>& hierarchy,
	int retr = RETR_LIST, int method = CHAIN_APPROX_SIMPLE, Point offset = Point(0, 0))
{
	/*
	*���ذ�����ͨ����Ķ�ֵͼ����
	*��ͨ������������������ܳ���
	*�������ټ����׿���MC
	*��ʹ���������ͨ����Ҳֻ���ص�һ��
	*/
	CvMat c_image = src;
	MemStorage storage(cvCreateMemStorage());
	CvSeq* _ccontours = 0;
	cvFindContours(&c_image, storage, &_ccontours, sizeof(CvContour), retr, method, CvPoint(offset));

	if (!_ccontours)
	{
		contours.clear();
		return -1;
	}
	Seq<CvSeq*> all_contours(cvTreeToNodeSeq(_ccontours, sizeof(CvSeq), storage));

	SeqIterator<CvSeq*> it = all_contours.begin();
	CvSeq* c = *it;
	return (int)c->total;
}

int CalculatePerimeterFast(BTNode bTNode, Mat srimg)
{
	/*
	*���ټ�����ͨ������ܳ�
	*/
	int height = srimg.rows;
	int width = srimg.cols;
	Mat m = Mat::zeros(height, width, CV_8UC1);

	//�����ͨ����
	for (int i =0; i<bTNode.pixelLocation.size(); i++)
		m.data[bTNode.pixelLocation[i]] = 255;

	vector<vector<Point>> contour;
	vector<Vec4i> hierarchy;
	return MyfindContours(m, contour, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_NONE, cv::Point());
}

double CalcuteCriterion_HeterogeneityChange(BTNode* t1, BTNode* t2, Mat srimg)
{
	/*
	*������׼��_3
	*�����Ա仯
	*������eCognition
	*/
	const double WS = 0.5;
	const double WFC = 0.25;
	const double WFS = 0.25;
	//�����ڵ�ϲ�����
	BTNode bTNode_temp;
	for (int i = 0; i<t1->pixelLocation.size(); i++)
		bTNode_temp.pixelLocation.push_back(t1->pixelLocation[i]);
	for (int i = 0; i<t2->pixelLocation.size(); i++)
		bTNode_temp.pixelLocation.push_back(t2->pixelLocation[i]);
	bTNode_temp.area = t1->area + t2->area;
	bTNode_temp.avgB = (t1->area*t1->avgB + t2->area*t2->avgB) / (t1->area + t2->area);
	bTNode_temp.avgG = (t1->area*t1->avgG + t2->area*t2->avgG) / (t1->area + t2->area);
	bTNode_temp.avgR = (t1->area*t1->avgR + t2->area*t2->avgR) / (t1->area + t2->area);
	
	//clock_t startTime,endTime; 
	//startTime = clock();
	//��һ���ڵ�
	double SD_1[3];
	if (t1->spectualStandDeviation_B != -1)
	{
		SD_1[0] = t1->spectualStandDeviation_B;
		SD_1[1] = t1->spectualStandDeviation_G;
		SD_1[2] = t1->spectualStandDeviation_R;
	}
	else
	{
		CalculateSpectralStandardDeviation(t1, SD_1, srimg);
		t1->spectualStandDeviation_B = SD_1[0];
		t1->spectualStandDeviation_G = SD_1[1];
		t1->spectualStandDeviation_R = SD_1[2];
	}
	//�ڶ����ڵ�
	double SD_2[3];
	if (t2->spectualStandDeviation_B != -1)
	{
		SD_2[0] = t2->spectualStandDeviation_B;
		SD_2[1] = t2->spectualStandDeviation_G;
		SD_2[2] = t2->spectualStandDeviation_R;
	}
	else
	{
		CalculateSpectralStandardDeviation(t2, SD_2, srimg);
		t2->spectualStandDeviation_B = SD_2[0];
		t2->spectualStandDeviation_G = SD_2[1];
		t2->spectualStandDeviation_R = SD_2[2];
	}
	//�ϲ��ڵ�
	double SD_12[3];
	CalculateSpectralStandardDeviation(&bTNode_temp, SD_12, srimg);

	//���ײ��������Ա仯
	double h_color = 0;
	for (int i = 0; i<3; i++)
		h_color += WS * abs(((t1->area+ t2->area)*SD_12[i] - t1->area*SD_1[i] - t2->area*SD_2[i]));
	//printf("h_color:%lf ", h_color);
	//endTime = clock();
	//printf("hcolor time:%lf, ", (double)(endTime - startTime) / CLOCKS_PER_SEC);


	//����bounding box�ܳ�
	//clock_t startTime_1,endTime_1;
	//startTime_1 = clock();
	int b[3] = {0};
	if (t1->boundingBoxPerimeter != -1)
		b[0] = t1->boundingBoxPerimeter;
	else
	{
		b[0] = CalculateBoundingBoxPerimete(t1, srimg);
		t1->boundingBoxPerimeter = b[0];
	}
	if (t2->boundingBoxPerimeter != -1)
		b[1] = t2->boundingBoxPerimeter;
	else
	{
		b[1] = CalculateBoundingBoxPerimete(t2, srimg);
		t1->boundingBoxPerimeter = b[1];
	}
	b[2] = CalculateBoundingBoxPerimete(&bTNode_temp, srimg);
	//endTime_1 = clock();
	//printf("bounding box time:%lf, ", (double)(endTime_1 - startTime_1) / CLOCKS_PER_SEC);

	//�����ܳ�
	//�����ڵ��ܳ�
	//clock_t startTime_2, endTime_2;
	//startTime_2 = clock();
	int l[3] = {0};
	if (t1->perimeter != -1)
		l[0] = t1->perimeter;
	else
	{
		l[0] = CalculatePerimeterFast(*t1, srimg);
		t1->perimeter = l[0];
	}
	if (t2->perimeter != -1)
		l[1] = t2->perimeter;
	else
	{
		l[1] = CalculatePerimeterFast(*t2, srimg);
		t2->perimeter = l[1];
	}
	l[2] = CalculatePerimeterFast(bTNode_temp, srimg);
	//endTime_2 = clock();
	//printf("perimeter time:%lf\n", (double)(endTime_2 - startTime_2) / CLOCKS_PER_SEC);
	//clock_t startTime_2, endTime_2;
	//startTime_2 = clock();
	//int width = srimg.cols;
	//int height = srimg.rows;
	//int* label = new int[width*height];
	//bool* visited = new bool[width*height];
	//for (int i = 0; i<height; i++)
	//	for (int j = 0; j<width; j++)
	//		label[i*width+j] = 0;
	//for (int i = 0; i<t1->pixelLocation.size(); i++)
	//	label[t1->pixelLocation[i]] = 1;
	//for (int i = 0; i<height; i++)
	//	for (int j = 0; j<width; j++)
	//		if (label[i*width+j] == 1)
	//			visited[i*width+j] = false;
	//		else
	//			visited[i*width+j] = true;	
	//int l[3] = {0};
	//if (t1->perimeter != -1)
	//	l[0] = t1->perimeter;
	//else
	//{
	//	CalcutePerimeter(t1->pixelLocation[0], label, visited, l[0], width, height);
	//	t1->perimeter = l[0];
	//}
	////����
	//for (int i = 0; i<height; i++)
	//	for (int j = 0; j<width; j++)
	//		label[i*width+j] = 0;
	//for (int i = 0; i<t2->pixelLocation.size(); i++)
	//	label[t2->pixelLocation[i]] = 1;
	//for (int i = 0; i<height; i++)
	//	for (int j = 0; j<width; j++)
	//		if (label[i*width+j] == 1)
	//			visited[i*width+j] = false;
	//		else
	//			visited[i*width+j] = true;
	//if (t2->perimeter != -1)
	//	l[1] = t2->perimeter;
	//else
	//{
	//	CalcutePerimeter(t2->pixelLocation[0], label, visited, l[1], width, height);
	//	t2->perimeter = l[1];
	//}
	////����
	//for (int i = 0; i<height; i++)
	//	for (int j = 0; j<width; j++)
	//		label[i*width+j] = 0;
	//for (int i = 0; i<t1->pixelLocation.size(); i++)
	//	label[t1->pixelLocation[i]] = 1;
	//for (int i = 0; i<t2->pixelLocation.size(); i++)
	//	label[t2->pixelLocation[i]] = 1;
	//for (int i = 0; i<height; i++)
	//	for (int j = 0; j<width; j++)
	//		if (label[i*width+j] == 1)
	//			visited[i*width+j] = false;
	//		else
	//			visited[i*width+j] = true;	
	//CalcutePerimeter(t1->pixelLocation[0], label, visited, l[2], width, height);
	//endTime_2 = clock();
	//printf("perimeter time:%lf\n", (double)(endTime_2 - startTime_2) / CLOCKS_PER_SEC);

	double h_smooth, h_compt;
	h_smooth = abs((t1->area + t2->area)*(double)l[2]/(double)b[2] - t1->area*(double)l[0]/(double)b[0] - t2->area*(double)l[1]/(double)b[1]);
	h_compt = abs((t1->area + t2->area)*(double)l[2]/sqrt((double)(t1->area + t2->area)) - t1->area*(double)l[0]/sqrt((double)t1->area) - t2->area*(double)l[1]/sqrt((double)t2->area));

	//������
	char outdatastr[] = "./data/hdata.txt";
	int flag = remove(outdatastr);
	FILE *fp;
	if((fp = fopen(outdatastr, "a")) == NULL)
	{
		printf("��hdata��¼�ļ�ʧ��\n");
		system("pause");
		exit(-1);
	}
	fprintf(fp, "sd: %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", SD_1[0], SD_1[1], SD_1[2], SD_2[0], SD_2[1], SD_2[2], SD_12[0], SD_12[1], SD_12[2]);
	fprintf(fp, "usual: %d %d %d %d %d %d %d %d %lf %lf %lf\n",t1->area, t2->area, l[0], l[1], l[2], b[0], b[1], b[2], h_color, h_smooth, h_compt);
	fclose(fp);

	//��ʱ�ͷ��ڴ�
	//delete label;
	//delete visited;

	return WS*h_color + WFS*h_smooth + WFC*h_compt;
}

double CalcuteCriterion_SpectualStandardDevation(BTNode* t1, BTNode* t2, Mat srimg)
{
	/*
	*������׼��_4
	*���ױ�׼��
	*/
	//�����ڵ�ϲ�����
	BTNode bTNode_temp;
	for (int i = 0; i<t1->pixelLocation.size(); i++)
		bTNode_temp.pixelLocation.push_back(t1->pixelLocation[i]);
	for (int i = 0; i<t2->pixelLocation.size(); i++)
		bTNode_temp.pixelLocation.push_back(t2->pixelLocation[i]);
	bTNode_temp.area = t1->area + t2->area;
	bTNode_temp.avgB = (t1->area*t1->avgB + t2->area*t2->avgB) / (t1->area + t2->area);
	bTNode_temp.avgG = (t1->area*t1->avgG + t2->area*t2->avgG) / (t1->area + t2->area);
	bTNode_temp.avgR = (t1->area*t1->avgR + t2->area*t2->avgR) / (t1->area + t2->area);
	//��һ���ڵ�
	double SD_1[3];
	if (t1->spectualStandDeviation_B != -1)
	{
		SD_1[0] = t1->spectualStandDeviation_B;
		SD_1[1] = t1->spectualStandDeviation_G;
		SD_1[2] = t1->spectualStandDeviation_R;
	}
	else
	{
		CalculateSpectralStandardDeviation(t1, SD_1, srimg);
		t1->spectualStandDeviation_B = SD_1[0];
		t1->spectualStandDeviation_G = SD_1[1];
		t1->spectualStandDeviation_R = SD_1[2];
	}
	//�ڶ����ڵ�
	double SD_2[3];
	if (t2->spectualStandDeviation_B != -1)
	{
		SD_2[0] = t2->spectualStandDeviation_B;
		SD_2[1] = t2->spectualStandDeviation_G;
		SD_2[2] = t2->spectualStandDeviation_R;
	}
	else
	{
		CalculateSpectralStandardDeviation(t2, SD_2, srimg);
		t2->spectualStandDeviation_B = SD_2[0];
		t2->spectualStandDeviation_G = SD_2[1];
		t2->spectualStandDeviation_R = SD_2[2];
	}
	//�ϲ��ڵ�
	double SD_12[3];
	CalculateSpectralStandardDeviation(&bTNode_temp, SD_12, srimg);
	//���ײ��������Ա仯
	double h_color = 0;
	for (int i = 0; i<3; i++)
		h_color += ((t1->area+ t2->area)*SD_12[i] - t1->area*SD_1[i] - t2->area*SD_2[i]);
	return h_color;
}


#endif