#include "lvf_scans.h"

extern int _FLAG_show_galaxy_;

/*
enum MENU_TYPE
{
        MENU_FRONT,
        MENU_SPOT,
        MENU_BACK,
        MENU_BACK_FRONT,
};
*/

//MENU_TYPE show = MENU_BACK_FRONT;

void menu(int item)
{
/*
        switch (item)
        {
        case MENU_FRONT:
        case MENU_SPOT:
        case MENU_BACK:
        case MENU_BACK_FRONT:
                {
                        show = (MENU_TYPE) item;
                }
                break;
        default:
                break;
        }
*/

    if(item) _FLAG_show_galaxy_=1;
    else _FLAG_show_galaxy_=0;

        glutPostRedisplay();

        return;
}


