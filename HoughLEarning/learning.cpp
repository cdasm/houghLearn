#include "learning.h"



auto learningMapForInstance(const vector<vector<int> >& votes,const vector<int>& positions)->vector<unordered_map<int,vector< pair<pair<int,int>,int>>>>
{
	assert(votes.size()>0);
	assert(votes[0].size()==rankNumber);

	//rank> position> category

	const int rankIndx=0;
	const int posIndx=1;
	const int categoryIndx=2;


	vector<vector<vector<vector<int> > > > allmaps(3);

	allmaps[rankIndx]=vector<vector<vector<int> > >(rankNumber,vector<vector<int>>(positionNumber,vector<int>(categoryNumber,0)));
	allmaps[posIndx]=vector<vector<vector<int> > > (positionNumber,vector<vector<int>>(categoryNumber,vector<int>(rankNumber,0))); 
	allmaps[categoryIndx]=vector<vector<vector<int> > > (categoryNumber,vector<vector<int>>(rankNumber,vector<int>(positionNumber,0))); 


	for (int i = 0; i < votes.size(); i++)
	{
		int pos=positions[i];
		for (int j = 0; j < votes[i].size(); j++)
		{
			int rank=j;
			int category=votes[i][j];
			++ allmaps[rankIndx] [rank][pos][category];
			++allmaps[posIndx][pos][category][rank];
			++allmaps[categoryIndx][category][rank][pos];
		}
	}

	vector<unordered_map<int,vector< pair<pair<int,int>,int>>>> result(3);

	for(int mi=0;mi<allmaps.size();++mi)
	for (int i = 0; i < allmaps[mi].size(); i++)
	{
		for (int j = 0; j < allmaps[mi][i].size(); j++)
		{
			for (int k = 0; k < allmaps[mi][i][j].size(); k++)
			{
				if (allmaps[mi][i][j][k]>0)
				{
					result[mi][i].push_back(make_pair(make_pair(j,k),allmaps[mi][i][j][k]));
				}
			}
		}
	}


	return result;

}




pair<MatrixXd,MatrixXd> funcAndJacobian(vector<unordered_map<int,vector< pair<pair<int,int>,int>>>>& rankposcatMap,const MatrixXd& weights,int trueCategory)
{
	//w=(a*x*x+b)/(c*x*x+d)
	//(2*x*(a*d - b*c))/(c*x^2 + d)^2
	int totalWeightNumber=rankNumber+positionNumber+categoryNumber-3;
	assert(weights.cols()==totalWeightNumber);

	const int rankIndx=0;
	const int posIndx=1;
	const int categoryIndx=2;

	vector<int> numbers(3);
	numbers[rankIndx]=rankNumber;
	numbers[posIndx]=positionNumber;
	numbers[categoryIndx]=categoryIndx;

	vector<int> startWeightIndx(3);
	startWeightIndx[0]=0-1;
	startWeightIndx[1]=rankNumber-1-1;
	startWeightIndx[2]=rankNumber-1+positionNumber-1-1;

	auto getWeight=[=](int _what,int _where)->double
	{
		if(_where==0)
			return 1;

		return weights(0,startWeightIndx[_what]+_where);
	};

	auto getWeightPos=[=](int _what,int _where)->int
	{
		assert(_where!=0);
			

		return startWeightIndx[_what]+_where;
	};

	//rank> position> category
	//category rank position
	//position category rank

	assert(rankposcatMap[categoryIndx].count(trueCategory));

	double goodWeight=0;
	double otherWeight=0;
	
	unordered_map<int,double> otherCoefficient;
	for (auto& cat:rankposcatMap[categoryIndx])
	{
		double curweight=0;
		for (int i = 0; i < cat.second.size(); i++)
		{
			double rankWeight=getWeight(rankIndx,cat.second[i].first.first);
			double posWeight=getWeight(posIndx,cat.second[i].first.second);
			double number=cat.second[i].second;
			curweight+=rankWeight*rankWeight*posWeight*posWeight*number;
		}
		otherCoefficient[cat.first]=curweight;
		double catWeight=getWeight(categoryIndx,cat.first);
		curweight*=catWeight*catWeight;

		if(cat.first==trueCategory)
			goodWeight+=curweight;
		else
			otherWeight+=curweight;
	}
	//w=(a*x*x+b)/(c*x*x+d)
	//(2*x*(a*d - b*c))/(c*x^2 + d)^2
	MatrixXd fvalue(1,1),jacvalue(MatrixXd::Zero( 1,totalWeightNumber));



	fvalue(0,0)=otherWeight/goodWeight;
	
	for (auto& cw:otherCoefficient)
	{
		if(trueCategory!=0)
		if(cw.first!=trueCategory)
		{
			int index=getWeightPos(categoryIndx,cw.first);
			jacvalue(0,index)= cw.second*2*getWeight(categoryIndx,cw.first)/goodWeight;
		}
		else
		{
			int index=getWeightPos(categoryIndx,cw.first);
			jacvalue(0,index)=-2*otherWeight/cw.second/pow(getWeight(categoryIndx,cw.first),3);
		}
	}


	//rank> position> category
	//category rank position
	//position category rank
		//w=(a*x*x+b)/(c*x*x+d)
	//(2*x*(a*d - b*c))/(c*x^2 + d)^2

	for (int i = 0; i < 2; i++)
	{
		for (auto& onemap:rankposcatMap[i])
		{
			if(onemap.first!=0)
			{
				double a(0.0),b(0.0),c(0.0),d(0.0);
				for (int j = 0; j < onemap.second.size(); j++)
				{
					if(i==rankIndx)
					{
						double posWeight=getWeight(posIndx,onemap.second[j].first.first);
						double categoryWeight=getWeight(categoryIndx,onemap.second[j].first.second);
						double number=onemap.second[j].second;
						//curweight+=rankWeight*rankWeight*posWeight*posWeight*number;

						double curweight=number*posWeight*posWeight*categoryWeight*categoryWeight;
						if(onemap.second[j].first.second==trueCategory)
						{
							c+=curweight;
						}
						else
						{
							a+=curweight;
						}

					}
					else if(i==posIndx)
					{
						double categoryWeight=getWeight(categoryIndx,onemap.second[j].first.first);
						double rankWeight=getWeight(rankIndx,onemap.second[j].first.second);
						double number=onemap.second[i].second;
						double curweight=rankWeight*rankWeight*categoryWeight*categoryWeight*number;
						if(onemap.second[j].first.first==trueCategory)
						{
							c+=curweight;
						}
						else
						{
							a+=curweight;
						}

					}

				}
				int index=getWeightPos(i,onemap.first);
				double x=getWeight(i,onemap.first);
				d=goodWeight-c*x;
				b=otherWeight-a*x;

				double t=(c*x*x + d);
				jacvalue(0,index)=(2*x*(a*d - b*c))/t/t;
			}
			
		}
	}
	

	return make_pair(fvalue,jacvalue);
}