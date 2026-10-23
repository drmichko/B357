#include "common.h"
int getinfoparms( int num, size_t *stab, size_t *orb, size_t *last )
{   char buffer[1024];
    FILE* srcFile = fopen( "info-3-5-7.txt", "r" );
    if ( ! srcFile ) {
	    perror( "info?" );
	    exit(1);
    }
    while ( fgets( buffer, 1024, srcFile ) ){
	    int tmp;
	    if ( sscanf( buffer, "#num=%d  stabSize=%ld orbits=%ld last=%ld", &tmp, stab, orb, last  ) == 4 ) {
		    if ( tmp == num ) {
		    	fclose( srcFile);
			    return 1;
		    }
	    }
    }
 fclose( srcFile);
 return 0;
}


boole getfullinfos( int num, size_t *stab, size_t *orb, size_t *last , aglGroup * grp )
{   char buffer[1024];
    FILE* srcFile = fopen( "info-3-5-7.txt", "r" );
    if ( ! srcFile ) {
            perror( "info?" );
            exit(1);
    }
    while ( fgets( buffer, 1024, srcFile ) ){
            int tmp;
            if ( sscanf( buffer, "#num=%d  stabSize=%ld orbits=%ld last=%ld", &tmp, stab, orb, last  ) == 4 ) {
                    if ( tmp == num ) {
			    size_t size;
			boole  h  = loadaglboolesize( srcFile, grp,  &size);
			assert( *stab ==  size );
                        fclose( srcFile);
                        return h;
                    }
            }
    }
 fclose( srcFile);
 return NULL;
}

uchar *loadZiptarget(char *dirname, int target, basis_t *b)
{
    FILE *srcFile;
    size_t lus;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;
    size = size / 8;
    uchar *res;
    sprintf(fn, "%s/%d-%d-%d-%d.zip", dirname, target, b->s, b->t, b->m);
    srcFile = fopen(fn, "r");
    if ( ! srcFile ) {
        perror( fn );
        exit( 1 );
    }
    res = calloc(size, sizeof(uchar));
    if (!res) {
        perror("allocation table");
        exit(1);
    }

    if (!srcFile) {
        perror(fn);
        exit(1);
    }
    lus = fread(res, 1, size, srcFile);
    if (lus != size) {
        printf("load problem:-(\n");
        printf("lus=%ld spaceSize=%ld\n", (long) lus, (long) size);
        exit(1);
    }
    fclose(srcFile);
    return res;
}




