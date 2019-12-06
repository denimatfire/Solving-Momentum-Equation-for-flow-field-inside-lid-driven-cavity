#include <iostream>
#include<math.h>
#include<fstream>
#include<time.h>

using namespace std;

int main()
{
    int i,j,n,m,iterr=0;
    float Re,b,x,y,l,h,sum=0,serror=10,werror=10,dx,dy,s1[101][101],w1[101][101],u[101][101],v[101][101],w0[101][101],s0[101][101],c1,c2,c3,c4;

    clock_t begin,end;
    double cpu_time_used;

    ofstream f1,f2,f3,f4;
    f1.open("VS_UV_DATA.plt");
    f3.open("y_VS_U.txt");
	f2.open("x_VS_V.txt");
	f4.open("SWerrorVSnoIteration.txt");

    cout<<"enter no of grid points along x axis : ";
    cin>>n;
    cout<<"\nenter no of grid points along y axis : ";
    cin>>m;
    cout<<"enter reynolds number :";
    cin>>Re;

    cout<<"\nenter the L and H :";
    cin>>l>>h;

    dx=l/(n-1);
    dy=h/(m-1);
    b=dx/dy;

    //u,v boundary conditions//
    for(i=0;i<n;i++){
        u[i][0]=0;
        v[i][0]=0;//bottom wall
        u[i][m-1]=1;
        v[i][m-1]=0;//top wall
    }

    for(j=0;j<m;j++){
        u[n-1][j]=0;//right wall
        v[m-1][j]=0;
        u[0][j]=0;//left wall
        v[0][j]=0;
    }
    //initial values of stream and vorticity
    for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			{
				s1[i][j]=0.01;
				s0[i][j]=0.01;
			}

	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			{
				w1[i][j]=0.5;
				w0[i][j]=0.5;
			}
			float cc=0.000001;
begin =clock();
do{
//solving stream fn and its error
    sum=0;
    for(i=1;i<n-1;i++){
        for(j=1;j<m-1;j++){
            s1[i][j]=(s0[i+1][j]+s0[i-1][j]+(b*b)*(s0[i][j+1]+s0[i][j-1])+w0[i][j]*dx*dx)/(2*(1+b*b));
            sum=sum+(s1[i][j]-s0[i][j]);
            s0[i][j]=s1[i][j];
        }
    }
    serror=sqrt((sum*sum)/((m-1)*(n-1)));
//update u v
for(i=1;i<n-1;i++){
    for(j=1;j<m-1;j++){
        u[i][j]=(s0[i][j+1]-s0[i][j-1])/(2*dy);
        v[i][j]=(s0[i-1][j]-s0[i+1][j])/(2*dx);
    }
}
//update vorticity BD
for(j=0;j<m;j++){
    w0[0][j]=(2*(s0[0][j]-s0[1][j]))/(dx*dx);//left wall
    w0[n-1][j]=(2*(s0[n-1][j]-s0[m-2][j]))/(dx*dx);//right wall
}
for(i=0;i<n;i++){
        w0[i][0]=(2*(s0[i][0]-s0[i][1]))/(dy*dy);//bottom wall
        w0[i][m-1]=(2*(s0[i][n-1]-s0[i][n-2]-1.0*dy))/(dy*dy);//top wall

}
//solve for vorticity and error
    sum=0;
for(i=1;i<n-1;i++){
    for(j=1;j<m-1;j++){
            c1=(1-(u[i][j]*dx*Re*0.5));
            c2=(1+(u[i][j]*dx*Re*0.5));
            c3=(b*b)*(1-(v[i][j]*dy*Re*0.5));
            c4=(b*b)*(1+(v[i][j]*dy*Re*0.5));
        w1[i][j]=(c1*w0[i+1][j]+c2*w0[i-1][j]+c3*w0[i][j+1]+c4*w0[i][j-1])/(2*(1+(b*b)));
        sum=sum+(w1[i][j]-w0[i][j]);
        w0[i][j]=w1[i][j];
    }
}
werror=sqrt((sum*sum)/((m-1)*(n-1)));
//updating the value of phi and omega
		for(i=0;i<n;i++)
			for(j=0;j<m;j++)
			{
				s0[i][j]=s1[i][j];
				w0[i][j]=w1[i][j];
			}


iterr++;
cout<<"iterration no- "<<iterr<<"\tserror- "<<serror<<"\twerror- "<<werror<<"\n";
f4<<"\n"<<iterr<<"\t"<<serror<<"\t"<<werror;
}while(serror>cc || werror>cc);

f1<<"ZONE"<<"\tI="<<n<<"\tJ="<<m;
    for(x=0,i=0;i<n;i++,x=x+dx){
        for(y=0,j=0;j<m;j++,y=y+dy){
            f1<<"\n"<<x<<"\t"<<y<<"\t"<<u[i][j]<<"\t"<<v[i][j]<<"\t"<<s0[i][j]<<"\t"<<w0[i][j];
        }
    }
for(y=0,i=0;i<n;i++,y=y+dy){
    f2<<"\n"<<y<<"\t"<<((v[i][51]+v[i][51])*0.5);
}
for(x=0,j=0;j<m;j++,x=x+dx){
    f3<<"\n"<<x<<"\t"<<((u[51][j]+u[50][j])*0.5);
}


end = clock();
cpu_time_used =((double)(end-begin))/CLOCKS_PER_SEC;
cout<<" \ntime used "<<cpu_time_used;


    return 0;
}
