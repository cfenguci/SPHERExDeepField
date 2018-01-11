#include <stdio.h>
#include <stdlib.h>
#include "math.h"

#include "string.h"

int main (int argc , char * argv [] )
{

        // argv[1] = 1.jpg

        FILE *file;
        char *buffer;
        unsigned long fileLen;

        //Open file
        file = fopen(argv[1], "rb");
        if (!file)
        {
                fprintf(stderr, "Unable to open file %s", argv[1]);
                return;
        }

        //Get file length
        fseek(file, 0, SEEK_END);
        fileLen=ftell(file);
        fseek(file, 0, SEEK_SET);

        printf("%e\n",fileLen);
        //Allocate memory
        buffer=(char *)malloc(fileLen+1);
        if (!buffer)
        {
                fprintf(stderr, "Memory error!");
                                fclose(file);
                return 1;
        }

       read(file,&buffer,sizeof(buffer));
       fclose(file);

       int i=0;

	while (i < sizeof(buffer))
	{
     		printf("%02X\n",(int)buffer[i]);
     		i++;
	}

        return 0;
}
