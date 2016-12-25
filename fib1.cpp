#include<iostream>
//#include<cstdlib>
//#include <tchar.h>
//#include <ncurses.h>
#include<cmath>
#include<iomanip>
#include<vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include<sstream>
#include <ctime>
#include <deque>
#include<time.h>
#include<stdio.h>
#include<map>
#include <bits/stdc++.h>
#include <utility>
#include<list>
#include<fstream>
#include<set>
#include<stdlib.h>
#define MIN -100009
//#define pb push_back
#define pp pair<int,int>
#define lp list<pair<int,int> >
using namespace std;
static int count=1;
//deckey

pair<int,int> p1;
namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}
std::stringstream ss;
using namespace std;

vector<string>v;
class Node{
	public:

	    int key;
		int rank;
        pair<int,int> kp;
		Node* parent;
		Node* child;
		Node* before;
		Node* after;

		bool mark;
		bool child_visit;

		Node(){};

		Node(pair <int,int> p1){
			child_visit = false;
            kp.first=p1.first;
            kp.second=p1.second;
            key=kp.second;
			rank=0;
			parent = child = NULL;
			before = after = NULL;
			mark = false;
		}

};

Node* dummy=NULL;
map<int,Node*>m;
//hash_map h<int,Node*>h;
class fib_heap{
	public:
		Node* min;
		int num_of_nodes;

	    fib_heap()
		{
			min=NULL;
			num_of_nodes=0;
		}

		Node* Insert( pair<int,int>  );
		Node* Delete_Min();
		Node* Meld(Node* x,Node* y);
		Node* fib_heap_link(Node* x,Node* y);
		Node* add_child(Node* x,Node* y);
		Node* cut(Node*z);
		Node* Decrease_key(int k,int v);
		void cascading_decrease_ranks(Node* y);
		Node* search(Node* x, int key);
		Node* search_aid(Node*x,int k);
		Node* search1(Node*x,int k);
		Node* search_auxillary(Node*x,int k);
		Node*Delete(int k);
};

Node* fib_heap::Insert( pair<int,int> p1)
{
	Node* x = new Node(p1);
	//key=x->kp.second;
  m[p1.first]=x;
  //  h.put(val,x);
   	if (min==NULL)
	  {
	  	min=x;
	  	num_of_nodes++;
	  	return x;
	  		//num_of_nodes++;

	  }

	num_of_nodes++;
	min=fib_heap_link(min,x);
	return min;
}

Node* fib_heap:: fib_heap_link(Node* x,Node* y)
{
	if(x->key > y->key)
	{
		Node*y_ptr=add_child(x,y);
		return y_ptr;
	}
	else
	{
        Node*x_ptr=add_child(y,x);
		return x_ptr;
    }
}

Node* fib_heap::add_child(Node* x,Node* y)
{
	Node*z=NULL;
	x->parent=y;
	z=y->child;
	x->before=NULL;
	x->after=z;
	if(z!=NULL)
	  z->before=x;
	y->child=x;
	return y;
}

Node* fib_heap::Meld(Node* min1,Node* min2)
{
	if(min1==NULL)
	  return min2;
	if(min2==NULL)
	  return min1;
	return fib_heap_link(min1,min2);

}



    Node* fib_heap:: Delete_Min()
    {
        Node*z=min;
        Node* y=NULL;
        Node* x=min->child;
    /*    Node*x_ptr=NULL;
    if(!z){return dummy;}
    if(!x){
            free(z);
    min=NULL;
            return dummy;
    }*/
     //   free(z);
           // z=NULL;
        int max_rank=0;
       // int num_of_nodes=100;
//int D = (int)log2(num_of_nodes) + 2;
      int D=10000;
        Node** A;
        A = (Node**)malloc(D*sizeof(Node*));
        for(int i = 0; i < D ; i++ )
            A[i] = NULL;
        while(x!=NULL)
        {
            y=x;
            x=x->after;
            Node*temp=NULL;
            Node*temp1=NULL;
         //   if(y==NULL){ break; return dummy;}
          //  while (y!=NULL && A[y->rank]!=NULL)
        while ( y!=NULL && A[y->rank]!=NULL)

            {


                if(y==NULL){
                    break;
                    }
                if(A[y->rank]->key > y->key){
                    temp=cut(A[y->rank]);
                    temp1=y;
                }
                else{

                //if(A[y->rank]->val > y->val){
                    temp=cut(y);
                    temp1=A[y->rank];
                }
                y=fib_heap_link(temp,temp1);
                A[y->rank]=NULL;
                y->rank=y->rank+1;
            }
            A[y->rank]=y;
            if(y->rank >max_rank)
            {
                max_rank=y->rank;
            }
        }

        for(int i = 0; i <= max_rank ; i++ )
        {
            if(A[i]!=NULL)
            {
                if(x==NULL)
                 x=A[i];
                else
                {
                  Node *temp,*temp1=NULL;
                  if(A[i]->key > x->key){
                    temp=cut(A[i]);
                    temp1=x;
                  }
                  else{

                //if(A[y->rank]->val > y->val){
                    temp=cut(x);
                    temp1=A[i];
                  }
                 //x=link(x,A[i]);
                 x=fib_heap_link(temp1,temp);
               //  A[i]=NULL;
                }
                 A[i]=NULL;;
            }

        }
            num_of_nodes--;
            min=x;
delete [] A;
           // free(A);
            //A=NULL;
            return min;
    }




Node* fib_heap:: cut(Node* x)
{
	Node*y=x->parent;
	Node *t=NULL;
	if(y->child==x)
	{
	  y->child=x->after;
    }

	if(x->before !=NULL)
	{
	  x->before->after=x->after;
          t=x->before;
	  x->before=NULL;
    }

	if(x->after !=NULL)
	{
	  x->after->before=t;
	  x->after=NULL;
    }
	return x;
}



Node* fib_heap::search_aid(Node*min,int k)
{
    Node*ptr=min;
    Node* z=NULL;
	z=search1(ptr,k);
	if(z!=NULL)
	 return z;
	if(!z)
	{
		Node*p1=ptr;
		while(p1 !=NULL)
		{
			if(p1->child_visit==true)
			{
			  p1->child_visit=false;
			  z=search_aid(p1->child,k);
			 if(z!=NULL)
			  return z;
		    }
			p1=p1->after;
		}
	}
	return dummy;
}

Node* fib_heap::search1(Node*z,int k)
{
	Node*ptr=z;
	while(ptr!=NULL){
	      if(ptr->child!=NULL)
	        ptr->child_visit=true;
       	 if(ptr->key==k)
		   return ptr;
		 ptr=ptr->after;
	    }
    return dummy;
}

Node* fib_heap ::search_auxillary(Node* z,int k)
{
	min->child_visit=true;
    if(k==z->key)
     return min;
    return search_aid(min->child,k);
}

deque<Node*> nQ;

Node* isThereQ(deque<Node*> *q,int n){
     for (long i=0; i<(long)q->size(); ++i)
     {
        if(q->at(i)->key==n)
          return q->at(i);
     }
     return NULL;
}


int plq=0;
void findLOrderQ(Node* r,int l,deque<Node*> *q,int &pl){
     if(r==NULL||isThereQ(q,r->key)!=NULL)
       return;

     q->push_back(r);
     Node *p=r;
     while(p->after!=NULL){
       findLOrderQ(p->after, l, q,pl);
       p=p->after;
     }
     if(r->child!=NULL){
       findLOrderQ(r->child, l+1,q,pl);
     }
}




Node* fib_heap::Decrease_key(int k,int v)
{
Node*x=m[k];
if(!x) return dummy;
	   x->kp.second=v;
	   x->key=v;
	   if(x==min)
	     return min;
	   min->mark=false;
	   cascading_decrease_ranks(x);
	   Node*x_ptr=cut(x);
       min=fib_heap_link(x_ptr,min);
       return min;
}

void fib_heap:: cascading_decrease_ranks(Node* y)
{
	y=y->parent;
	if(y!=NULL)
	{
	if(y->rank > 0)
	 y->rank=y->rank-1;
	y->mark=!y->mark;
	if(y->mark==false)
	 cascading_decrease_ranks(y);
    }
    return;
}

Node* fib_heap::Delete(int k)
{
	Node*z=Decrease_key(k,MIN);
	if(z==min && z->child==NULL){
		//free(z->child);
		//z->child=NULL;
	 	min=NULL;
	 	return dummy;
	 }
	if(!z)
	 {
	 	return dummy;
	 }

	 Node*min=Delete_Min();
 	 return min;
}


class NodeLOrder{
     public:
     Node *n;
     pair<int,int> l;
     NodeLOrder(Node *n,pair<int,int> l){
         this->n = n;
         this->l = l;
     }
};

bool isThere(deque<NodeLOrder*> *q,Node *n){
     for (long i=0; i<(long)q->size(); ++i)
     {
        if(q->at(i)->n==n)
          return true;
     }
     return false;
}


int pl=0;
void findLOrder(Node* r,int l,deque<NodeLOrder*> *q,int &pl){
     if(r==NULL||isThere(q,r))
       return;

     q->push_back(new NodeLOrder(r,make_pair(l,++pl)));
     Node *p=r;
     while(p->after!=NULL){
       findLOrder(p->after, l, q,pl);
       p=p->after;
     }
     if(r->child!=NULL){
       findLOrder(r->child, l+1,q,pl);
     }
}

NodeLOrder* isThereN(deque<NodeLOrder*> *q,Node *n){
     for (long i=0; i<(long)q->size(); ++i)
     {
        if(q->at(i)->n==n)
          return q->at(i);
     }
     return NULL;
}

deque<NodeLOrder*> names;




void traverseLevelOrder(Node *r){
         cout<<endl<<endl;
         if(!r)return;
         cout<<r->key<<"("<<r->rank<<")"<<","<<endl;

         Node *p=r->child;
         if(p==NULL)
           return;
         deque<Node*> q;
         q.push_back(p);

         while(p->after!=NULL){
           q.push_back(p->after);
           p=p->after;
         }
         int l=2;
         while(!q.empty()){
           Node *p=q.front();
           if(p->child!=NULL){
             p=p->child;
             q.push_back(p);
             while(p->after!=NULL){
               q.push_back(p->after);
               p=p->after;
             }
           }
           NodeLOrder *t=isThereN(&names,q.front());
           if(t!=NULL && t->l.first!=l){
               l=t->l.first;
               cout<<endl;
           }
           Node *b=NULL,*a=NULL;
            if(q.front()->before!=NULL)
              b=q.front()->before;
            if(q.front()->after!=NULL)
              a=q.front()->after;
                  if(b!=NULL)
                    cout<<b->key;
                  else
                    cout<<'N';
                  cout<<"<--*"<<q.front()->key<<"("<<q.front()->parent->key<<")"<<"("<<q.front()->rank<<")"<<"*-->";
                  if(a!=NULL)
                    cout<<a->key<<",";
                  else
                    cout<<'N'<<",";
                //}
                q.pop_front();
           //}

         }
         cout<<endl;
}

//int d[400000] = {};
int main()
{
	fib_heap* f1 = new fib_heap();
    int n,c;
	 ifstream infile;
	infile.open("/home/murukessan/Documents/Venkat/new/US-d.BAY.dat");
 // infile.open("C:\\Users\\kumar\\Downloads\\small.dat");

   if (!infile)
      exit ( EXIT_FAILURE );

   int m,w,v,v1,w1;

    infile >> n >> m;
    int *d;
   //d = (int*)malloc (sizeof(int) * n);
d=new int[n];
  // arr = (int*)malloc (sizeof(int) * n);

  //  d[0]=9999;
 // int* d = new int[1000000];
d[0]=-129999;
//arr[0]=-8922;
     for(int l=1;l<=n+1;l++)
    {
        d[l] =  9999999;
   //    arr[l]=-8922;
    }

   // vector< pair<int, int>  > adjacencyList[n + 1];

   vector< list< pair<int, int> >  >  adjacencyList(n+1);
//adjacencyList.reserve(n)
                         // for(int p=1;p<=n;++p) adjacencyList[p].reserve(1000);
fprintf(stdout,"n %d, m %d\n", n, m);

d[1]=0;
//Node*q=f1->Insert( make_pair(1,d[1]));
std::set<int>s1;
set <int> s2;
set<int> result;
for ( int k = 1; k <= m; k++ ) {
infile >> v >> w >> c;
adjacencyList[v].push_back(make_pair(w, c));
//arr[v]=v;
s1.insert(v);
s2.insert(w);
//Node*q1=f1->Insert(make_pair(v,d[v]));
}
infile.close();

std::set_difference(s2.begin(), s2.end(), s1.begin(), s1.end(),
    std::inserter(result, result.end()));

set<int>::iterator it1;
  for (it1 = s1.begin(); it1 != s1.end(); it1++) {
           Node*q11=f1->Insert(make_pair(*it1,d[*it1]));
  }


set<int>::iterator it;
  for (it = result.begin(); it != result.end(); it++) {
        cout<<*it<<endl;
           Node*q11=f1->Insert(make_pair(*it,d[*it]));
  }


  //	infile1.open("C:\\Users\\kumar\\Downloads\\ex.dat");

/*
  	  	infile1.open("C:\\Users\\kumar\\Downloads\\small.dat");

  	infile1>>n>>m;
//  infile1.open();
  //for (it = result.begin(); it != result.end(); it++){
  for ( int k1 = 1; k1 <= m; k1++ ){
    infile1 >>v >>w >> c;
    if (result.count(w)) {
            cout<<w<<endl;
 adjacencyList[v].push_back(make_pair(w, c));
 // adjacencyList[w].push_back();

    }}
infile1.close();

*/



/*
for(int t=1;t<=m;++t){
        lp::iterator it=adjacencyList[t].begin();
for(;it!=adjacencyList[t].end();++it)
    {cout<<t<<" "<<(*it).first<<" "<<(*it).second<<endl;
}}
*/


for (int t = 1; t < adjacencyList.size(); ++t) {
        printf("adjacencyList[%d] ", t);

        list< pair<int, int> >::iterator itr = adjacencyList[t].begin();

        while (itr != adjacencyList[t].end()) {
            printf(" -> %d(%d)", (*itr).first, (*itr).second);
            ++itr;
        }
        printf("\n");
    }


/*
 int pl=0;
                                        names.clear();
                                        findLOrder(f1->min,1,&names,pl);
                                        for (long i=0; i<(long)names.size(); ++i)
                                        {
                                          cout << names.at(i)->n->key << "," << names.at(i)->l.first << " ," << names.at(i)->l.second << endl;
                                        }

					traverseLevelOrder(f1->min);*/

while(f1->min != NULL){
      int u=f1->min->kp.first;
//cout<<f1->num_of_nodes<<endl;
      if(f1->min->child==NULL)
					{
					    f1->min=NULL;
                //        free(f1);
						fib_heap* f1 = new fib_heap();
						cout<<"Heap is empty. So no nodes to delete1\n"<<endl;
						cout<<"\n-----------------------------------------------------------------------------------------------------------\n";
					    cout<<"-----------------------------------------------------------------------------------------------------------\n";
					    return 0 ;
					}
					if(f1->num_of_nodes==1)
					{
					    f1->min=NULL;
                        //free(f1);
						fib_heap* f1 = new fib_heap();
						cout<<"Heap is empty. So no nodes to delete\n"<<endl;
						cout<<"\n-----------------------------------------------------------------------------------------------------------\n";
					    cout<<"-----------------------------------------------------------------------------------------------------------\n";
					    break;
					}
					else{
					Node* z1=f1->Delete_Min();
					}

  /*    pl=0;
                                        names.clear();
                                        findLOrder(f1->min,1,&names,pl);
                                        for (long i=0; i<(long)names.size(); ++i)
                                        {
                                          cout << names.at(i)->n->key << "," << names.at(i)->l.first << " ," << names.at(i)->l.second << endl;
                                        }*/
//if(f1->min) break;


//traverseLevelOrder(f1->min); //******uncomment

     // int size = adjacencyList[u].size();
    //  for(int i1 = 0 ; i1 < size ; i1++)
      //  {
            lp::iterator it3=adjacencyList[u].begin();
            for(;it3!=adjacencyList[u].end();++it3)
{
            v1 = (*it3).first;
            w1 = (*it3).second;


                if(d[v1] > d[u] + w1)
{
              //  int y=d[v1];
                 d[v1]=d[u] + w1;
                 Node*t=f1->Decrease_key(v1,d[v1]);

              //   traverseLevelOrder(t); //******uncomment
                }
}
//traverseLevelOrder(f1->min);
     // }
}

// for(int i2=1; i2<=m; i2++) printf("Node %d, min weight = %d\n", i2, d[i2]);

ofstream out_data("/home/murukessan/Documents/Venkat/new/filename.dat");
for(int i2=1; i2<n; i2++)
{
    out_data<<"node "<<i2<<" with min-weight "<<d[i2]<<endl;
}
delete [] d;
out_data.close();
/*     for(int i=0;i<n;i++){
    if(f1->min == NULL)
					{
						cout<<"Heap is empty. So no nodes to delete\n";
						cout<<"\n-----------------------------------------------------------------------------------------------------------\n";
     					cout<<"-----------------------------------------------------------------------------------------------------------\n";
						break;
					}

					if(f1->num_of_nodes==1)
					{
					    f1->min=NULL;
                        //free(f1);
						fib_heap* f1 = new fib_heap();
						cout<<"Heap is empty. So no nodes to delete1\n"<<endl;
						cout<<"\n-----------------------------------------------------------------------------------------------------------\n";
					    cout<<"-----------------------------------------------------------------------------------------------------------\n";
					    break;
					}
					else
					{
					Node* z=f1->Delete_Min();
				 	if(!f1->min){
				        	//cout<<"empty"<<endl<<endl;
				        	break;
						}}}

*/
				     return 0;

        }


