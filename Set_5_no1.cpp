
 /* *******************************************************
 * This program reads a matrix, prints it, computes       *
 * determinant and transpose of the matrix and finds the  *
 * sum, difference and product of matrices using pointers.*
 * Written by Ishaan Arora , 6-22 Sep 2018                *
 **********************************************************/

 /**Currently this program can perform operations on at most two matrices,
      but the number of matrices can be extended using classes.**/

#include<iostream>
#include<malloc.h>
#include<math.h>
using namespace std;

double** memory(double** pointer, int frow,int fcol);
double** identity(int frow, int fcol);
double** matread( double**ptr,int frow, int fcol);
void matprint(double**pt,int frow, int fcol,int x);
double** transpose(double** ptr,int frow,int fcol);
int nonzero(double **ptr,int frow,int i,int j);                 // used for  row pivoting
double** swap(double**detptr,int si,int sj,int frow,int fcol);  // used for row pivoting
double det(double** ptr,int frow,int fcol);
double** matadd(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1);
double** matsub(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1);
double** matmul(double **ptr1,double **ptr2,int frow,int fcol,int frow1,int fcol1);
double** option1(double**ptr,int frow,int fcol);
double** option2(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1);

double** memory(double** pointer, int frow,int fcol)
{
    pointer=(double**)malloc(sizeof(double)*frow);
    for(int i=0;i<frow;i++)
    {
        *(pointer+i)=(double*)malloc(sizeof(double)*fcol);
    }
    return pointer;
}

double** identity(int frow, int fcol)
{   double **id;
    id=memory(id,frow,fcol);
    for(int i=0;i<frow;i++)
    {
        for(int j=0;j<fcol;j++)
        {
            if(i==j)
                *(*(id+i)+j)=1;
            else
                *(*(id+i)+j)=0;
        }
    }
    return id;
}

double** matread( double**ptr,int frow, int fcol)
{
    ptr=memory(ptr,frow,fcol);

    cout<<"\nNow enter the elements : ";
    for(int i=0;i<frow;i++)
    {
        for(int j=0;j<fcol;j++)
        {
            cin>>*(*(ptr+i)+j);
        }
    }
    return ptr;         //this return is the most important part of this program.
}


void matprint(double**pt,int frow, int fcol,int x)
{
    cout<<"\n";
    for(int i=0;i<frow;i++)
    {
        for(int j=0;j<fcol;j++)
        {
            cout<<*(*(pt+i)+j)<<"\t";
        }
        cout<<"\n";
    }
    if(x!=0)
    {
    cout<<"\nIs this the matrix you entered ?\n[Y/N] :";
    char c;
    cin>>c;
    if(c!='y'&&c!='Y')
        exit(1);
    }
}

double** transpose(double** ptr,int frow,int fcol)
{
    double** tpose;
    tpose=memory(tpose,fcol,frow);

    for(int i=0;i<fcol;i++)
    {
        for(int j=0;j<frow;j++)
        {
            *(*(tpose+i)+j)=*(*(ptr+j)+i);
        }
    }
    matprint(tpose,fcol,frow,0);
    return tpose;
}

int nonzero(double **ptr,int frow,int i,int j)
{
    while(*(*(ptr+i)+j)==0)
    {
        i++;
        if(i==frow-1)
            break;

    }
    //cout<<i;
    return i;
}

double** swap(double**detptr,int si,int sj,int frow,int fcol)
{
            for(int j=0;j<fcol;j++,sj++)
                {
                    *(*(detptr+0)+j)=*(*(detptr+0)+j)+*(*(detptr+si)+sj);
                    *(*(detptr+si)+sj)=*(*(detptr+0)+j)-*(*(detptr+si)+sj);
                    *(*(detptr+0)+j)=*(*(detptr+0)+j)-*(*(detptr+si)+sj);
                     //matprint(detptr,frow,fcol);
                }

    return detptr;
}

double det(double** ptr,int frow,int fcol)
{
    double dete=1;
    double **detptr;
    detptr=memory(detptr,frow,fcol);

    for(int i=0;i<frow;i++)
        {
            for(int j=0;j<fcol;j++)
            {
               *(*(detptr+i)+j)=*(*(ptr+i)+j);
            }
        }

    if(frow==fcol)
    {
        if(*(*(detptr+0)+0)==0)
        {
            int swapi;
            swapi=nonzero(detptr,frow,0,0);
            detptr=swap(detptr,swapi,0,frow,fcol);
            dete*=-1;
        }


     for(int k=0;k<fcol;k++)
        {
                    int count=0,countcol=0;
                    for(int b=0;b<frow;b++)
                        {
                            for(int a=0;a<fcol;a++)
                            {
                            if(*(*(detptr+b)+a)==0)
                                count++;
                            if(*(*(detptr+a)+b)==0)
                                countcol++;

                            }
                        }
                    if(count==fcol||countcol==frow)
                    {
                        dete=0;
                        break;
                    }
            for(int i=frow-1;i>k;i--)
            {
                double ratio=double((*(*(detptr+i)+k)))/double((*(*(detptr+k)+k)));

                for(int j=0;j<fcol;j++)
                {
                    //cout<<"\nRatio:"<<ratio;
                    //cout<<"\nTo sub :"<<ratio*(*(*(detptr+k)+j));
                    //cout<<"\nk:"<<k<<'\t'<<"i:"<<i<<'\t'<<"j:"<<j;
                    *(*(detptr+i)+j)= *(*(detptr+i)+j)-(ratio)*(*(*(detptr+k)+j));

              // matprint(detptr,frow,fcol,0);
              // remove these upper comments to see the process of computing determinant

                }

                int count=0, x=1;		// checks whether the matrix has become upper triangular or not
                    for(int l=0;l<(frow-1);l++ )
                        {
                            for(int m=x;m<fcol;m++)
                                {
                                    if(*(*(detptr+m)+l)!=0)
                                        count=1;
                                }
                                        x++;
                        }
                            if(count==0)
                                break;

                    count=0;			// checks whether all elements in a row are 0
                    int countcol=0;
                    for(int b=0;b<frow;b++)
                        {
                            for(int a=0;a<fcol;a++)
                            {
                            if(*(*(detptr+b)+a)==0)
                                count++;
                            if(*(*(detptr+a)+b)==0)
                                countcol++;

                            }
                        }
                    if(count==fcol||countcol==frow)
                    {
                        dete=0;
                        break;
                    }
            }

        }

        for(int i=0;i<frow;i++)
        {
            for(int j=0;j<fcol;j++)
            {
                if(i==j)
                    dete=dete*(*(*(detptr+i)+j));
            }
        }
        cout<<dete;
        return dete;

    }
}
/*						// this part of code is still in development phase
double** inverse (double** ptr,int frow,int fcol)
{
    double**uptr;
    uptr=memory(uptr,frow,fcol);

    for(int i=0;i<frow;i++)
        {
            for(int j=0;j<fcol;j++)
            {
               *(*(uptr+i)+j)=*(*(ptr+i)+j);
            }
        }

        double**lptr;
        lptr=memory(lptr,frow,fcol);

    for(int i=0;i<frow;i++)
        {
            for(int j=0;j<fcol;j++)
            {
               *(*(lptr+i)+j)=*(*(ptr+i)+j);
            }
        }

//if(det(ptr,frow,fcol)!=0)

        if(*(*(uptr+0)+0)==0)
        {
            int swapi;
            swapi=nonzero(uptr,frow,0,0);
            uptr=swap(uptr,swapi,0,frow,fcol);

        }

     for(int k=0;k<fcol;k++)                                /* produces an upper triangular matrix */
  /*      {
            for(int i=frow-1;i>k;i--)
            {
                double ratio=double((*(*(uptr+i)+k)))/double((*(*(uptr+k)+k)));
                //cout<<"\nRatio :"<<double((*(*(detptr+i)+k)))/double((*(*(detptr+(i-1))+k)));
                for(int j=0;j<fcol;j++)
                {
                   // cout<<"\nRatio:"<<ratio;
                   // cout<<"\nTo sub :"<<ratio*(*(*(inptr+k)+j));
                   // cout<<"\nk:"<<k<<'\t'<<"i:"<<i<<'\t'<<"j:"<<j;
                    *(*(uptr+i)+j)= *(*(uptr+i)+j)-(ratio)*(*(*(uptr+k)+j));
                }

                int count=0, x=1;
                    for(int l=0;l<(frow-1);l++ )
                        {
                            for(int m=x;m<fcol;m++)
                                {
                                    if(*(*(uptr+m)+l)!=0)
                                        count=1;
                                }
                                        x++;
                        }
                            if(count==0)
                                break;
            }
        }
                matprint(uptr,frow,fcol,0);


                                    /** produces a lower triangular matrix **/
     /*    lptr=identity(frow,fcol);
         for(int i=frow-1;i>0;i--)
         {
             *(*(lptr+i)+0)=(*(*(ptr+i)+0))/(*(*(uptr+0)+0));
         }


         matprint(lptr,frow,fcol,0);
        matprint(matmul(lptr,uptr,frow,fcol,frow,fcol),frow,fcol,0);
        double** id;
        id=identity(frow,fcol);
        matprint(id,frow,fcol,0);
}*/

double** matadd(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1)
{
    if(frow1==frow&&fcol==fcol1)
    {
        double**sumptr;
        sumptr=memory(sumptr,frow,fcol);

        for(int i=0;i<frow;i++)
        {
            for(int j=0;j<fcol;j++)
            {
                *(*(sumptr+i)+j)= *(*(ptr1+i)+j) + *(*(ptr2+i)+j);
            }
        }
        matprint(sumptr,frow,fcol,0);
        return sumptr;
    }

    else
        cout<<"\n(o_o) Sorry, Matrix 1 and 2 can't be added (>_<).\n";
}

double** matsub(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1)
{
     if(frow1==frow&&fcol==fcol1)
    {
        double** subptr;
        subptr=memory(subptr,frow,fcol);

        char op;
        cout<<"\nDo you want to subtract :\na ;Matrix 1 from Matrix 2\tor\nb :Matrix 2 from Matrix 1\n";
        cout<<"Enter your choice : ";
        cin>>op;
        if(op=='a'||op=='A')
        {
            for(int i=0;i<frow;i++)
            {
                for(int j=0;j<fcol;j++)
                {
                    *(*(ptr1+i)+j)=-*(*(ptr1+i)+j);
                    *(*(subptr+i)+j)=*(*(ptr1+i)+j)+*(*(ptr2+i)+j);
                }
            }
        }

        else
        {   if(op!='b'&&op!='B')
            cout<<"I'll take that as a 'b' . :P";

            for(int i=0;i<frow;i++)
            {
                for(int j=0;j<fcol;j++)
                {
                    *(*(ptr2+i)+j)=-*(*(ptr2+i)+j);
                    *(*(subptr+i)+j)=*(*(ptr1+i)+j)+*(*(ptr2+i)+j);
                }
            }
        }

        matprint(subptr,frow,fcol,0);
        return subptr;
    }
    else
        cout<<"\n\(o_o)/ Sorry, Matrix 1 and 2 can't be subtracted (;-;).\n";
}

double ** matmul(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1)
{
    double**mul;
    char op;
    cout<<"\nDo you want to multiply :\na :Matrix 1 to Matrix 2 \t or\nb :Matrix 2 to Matrix 1\n";
    cout<<"Enter your choice :";
    cin>>op;
    if(op=='a'||op=='A')
    {   if(fcol==frow1)
        {
            mul=memory(mul,frow,fcol1);

            for(int i=0;i<frow;i++)
            {
                for(int k=0;k<fcol1;k++)
                {
                        int sum=0;
                    for(int j=0;j<fcol;j++)
                    {
                        sum=sum+(*(*(ptr1+i)+j))*(*(*(ptr2+j)+k));
                    }
                        *(*(mul+i)+k)=sum;
                }
            }
            matprint(mul,frow,fcol1,0);
        }
    else
        cout<<"\n(>_<) Sorry, Matrix 1 and 2 can't be multiplied (;-;).\n";
    }

    else if(op=='b'||op=='B')
    {   if(fcol1==frow)
        {
            mul=memory(mul,frow1,fcol);
            for(int i=0;i<frow1;i++)
            {
                for(int k=0;k<fcol;k++)
                {
                        int sum=0;
                    for(int j=0;j<fcol1;j++)
                    {
                        sum=sum+(*(*(ptr1+i)+j))+(*(*(ptr2+j)+k));
                    }
                        *(*(mul+i)+k)=sum;
                }
            }
            matprint(mul,frow1,fcol,0);
        }
        else
            cout<<"(>_<) Sorry, Matrix 2 and 1 can't be multiplied (;-;).";
    }
return mul;
}

double** option1(double**ptr,int frow,int fcol)
{
    double** result;
    cout<<"\nSelect an operation on matrix :\na :Transpose the matrix \t b :Find determinant of matrix\n";
    cout<<"Enter your choice :";
    char c;
    cin>>c;
    cout<<"\nThe result is :\n";
    if(c=='a'||c=='A')
    {
        result=transpose(ptr,frow,fcol);
        return result;
    }
    else if(c=='b'||c=='B')
    {
        det(ptr,frow,fcol);
       // return result;
    }
    /*else if(c=='c'||c=='C')
    {
        inverse(ptr,frow,fcol);
       // return result;
    }*/
}

double** option2(double**ptr1,double**ptr2,int frow,int fcol,int frow1,int fcol1)
{       double** result;
    cout<<"\nSelect an operation on matrices :\na :Add two matrices\tb :Subtract two matrices"<<
    "\nc :Multiply two matrices\n";
    cout<<"Enter your choice :";
    char c;
    cin>>c;
    cout<<"\nThe result is :\n";
    if(c=='a'||c=='A')
    {
        result=matadd(ptr1,ptr2,frow,fcol,frow1,fcol1);
        return result;
    }
    else if(c=='b'||c=='B')
    {
        result=matsub(ptr1,ptr2,frow,fcol,frow1,fcol1);
        return result;
    }
    else if(c=='c'||c=='C')
    {
        result=matmul(ptr1,ptr2,frow,fcol,frow1,fcol1);
        return result;
    }
}

int main()
{
    int no;
    cout<<"(Note : This program is capable of handling upto only two matrices)\n" ;
    cout<<"\nYou want to perform operation(s) on how many no. of matrices?\nEnter the number :";
    cin>>no;
    double**p,**result,row=0,col=0;
    if(no==1)
    {
        cout<<"\nHi!";
        cout<<"\nEnter the no of rows :";
        cin>>row;
        cout<<"Enter the no of columns :";
        cin>>col;
        p=matread(p,row,col);
        matprint(p,row,col,1);
        result=option1(p,row,col);
    }

    if(no==2)
    {
        cout<<"\nHi!";
        cout<<"\nEnter the no of rows :";
        cin>>row;
        cout<<"Enter the no of columns :";
        cin>>col;
        p=matread(p,row,col);
        matprint(p,row,col,1);
        int row2=0,col2=0;
        cout<<"\n\nFor second matrix\n";
        cout<<"\nEnter the no of rows :";
        cin>>row2;
        cout<<"Enter the no of columns :";
        cin>>col2;

        double**p1;
        p1=matread(p1,row2,col2);
        matprint(p1,row2,col2,1);
        result=option2(p,p1,row,col,row2,col2);
    }
}
