#include "raj_code.h"
#include "lvf_scans.h"

struct_file_dims * get_file_dims(string fileName)
{
  FILE *dataFile;
  char customer[1024];
  int nRow = 0, nCol = 0, n=0;
  char *line;
  char entry[1024]; 
 
  struct_file_dims *p=new struct_file_dims;

  p->filename=fileName;


  if((dataFile = fopen(fileName.c_str(), "r")) == NULL) 
  {
    printf("Error Opening File.\n");
    exit(1);
  }

  while((line = fgets(customer, sizeof(customer), dataFile)) != NULL )
  {
    while(sscanf(line, "%1024s%n", entry, &n ) == 1 )
    {
      nCol++;
      line += n;
    }

    nRow++;
    //printf("%d %d\n",nRow, nCol);
    //printf("%s\n",customer);

    p->nrow=nRow;
    p->ncol=nCol;

    nCol = 0;
  }
 
  fclose(dataFile);



  return p;
}

struct_data_table * load_data_table(string file)
{
  int i,j,k;
  struct_file_dims * fd=get_file_dims(file);
  int nx,ny;

  nx=fd->nrow;
  ny=fd->ncol;

  struct_data_table * dt=new struct_data_table;
  dt->nrow=nx;
  dt->ncol=ny;
  dt->filename=fd->filename;
  dt->table=create_2d_grid<double>(nx,ny);


  double *col=new double [ny];
  ifstream in;
  in.open(file.c_str());
  for(i=0;i<nx;i++) 
  {
    for(j=0;j<ny;j++) in>>col[j];
    for(j=0;j<ny;j++) dt->table[i][j]=col[j];
  }
  in.close();

  cout<<"Loaded "<<file<<endl;
  delete col;
  return dt;

}


void get_grid(string path, string name)
{

  struct_data_table *dt=load_data_table(path+name);
  int nrow=dt->nrow;
  cout<<"There are "<<nrow<<" entries"<<endl;
  int i,j;
  string file;


  double xmin=-FLAT_PATCH,xmax=FLAT_PATCH;
  double ymin=-FLAT_PATCH,ymax=FLAT_PATCH;
  double dx=FLAT_SPACE,dy=FLAT_SPACE;

  int numx=(xmax-xmin)/dx+1;
  int numy=(ymax-ymin)/dy+1;
  double x,y;

  double **pole_count=create_2d_grid<double>(numx,numy);


  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    pole_count[i][j]=0.0;
  }

  double theta,phi;
  for(i=0;i<nrow;i++)
  {
    double xx=dt->table[i][0];
    double yy=dt->table[i][1];
    double nhits=dt->table[i][2];
    int idx=(xx-xmin)/dx;
    int idy=(yy-ymin)/dy;
    pole_count[idx][idy]+=nhits;
   }


  ofstream out;
  //string null="/home/cfeng/computing/spherex/LVFscans/output/framedata/daily/";

  out.open((path+name+"-grid").c_str());
  for(i=0;i<numx;i++) for(j=0;j<numy;j++) 
  {
    out<<i<<"\t"<<j<<"\t"<<pole_count[i][j]<<endl;
    //cout<<i<<"\t"<<j<<"\t"<<pole_count[i][j]<<endl;
  } 
  out.close();

  delete_2d_grid<double>(pole_count,numx,numy);
  delete dt;

}
