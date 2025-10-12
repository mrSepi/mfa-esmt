#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <glpk.h>

typedef struct {
    double x;
    double y;
} Point;
Point *newPoint(double x, double y) {
    Point *res = malloc(sizeof (Point));
    res->x = x;
    res->y = y;
    return res;
}

typedef struct {
    int a;
    int b; /* indici dei due estremi */
    double len;
} Edge;
Edge *newEdge(int a, int b, double len) {
    Edge *res = malloc(sizeof (Edge));
    res->a = a;
    res->b = b;
    res->len = len;
    return res;
}

typedef struct node Node; /* per union-find */
struct node {
    int i; /* indice del vertice */
    Node *father;
    int height;
};
Node *newNode(int i) {
    Node *res = malloc(sizeof (Node));
    res->i = i;
    res->father = NULL;
    res->height = 0;
}

typedef struct item SteinerItem; /* per steinerization finale (evitare di ricalcolare mst) */
struct item {
    int index;
    Point *sp;
    SteinerItem *next;
};
SteinerItem *newItem(int index, Point *sp, SteinerItem *next) {
    SteinerItem *res = malloc(sizeof (SteinerItem));
    res->index = index;
    res->sp = sp;
    res->next = next;
    return res;
}

typedef struct elem Elem;
struct elem {
    Point *p; /* point */
    double w; /* weight */
    Elem *next;
};
Elem *newElem(Point *p, double w, Elem *next) {
    Elem *res = malloc(sizeof (Elem));
    res->p = p;
    res->w = w;
    res->next = next;
    return res;
}

double dist(Point *a, Point *b) {
    return sqrt(pow(a->x-b->x, 2) + pow(a->y-b->y, 2));
}


int compareEdgesByLen(const void *elem1, const void *elem2) {
    Edge *edge1 = *(Edge**) elem1; /* gli argomenti sono puntatori agli elementi dell'array */
    Edge *edge2 = *(Edge**) elem2;
    if (edge1->len < edge2->len) return -1;
    if (edge1->len > edge2->len) return  1;
    return 0;
} /* to sort edges for emst kruskal */
void emst(Point *points[], int n, bool adjacency[]) {
    /* adjacency matrice nxn */
    for (int i=0; i<n*n; i++)
        adjacency[i] = false; /* reset matrix: O(n^2)*/
    int numEdges = n*(n-1)/2;
    Edge **edges = malloc(sizeof (Edge*) * numEdges);
    int edgeIndex = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (j > i) edges[edgeIndex++] = newEdge(i,j,dist(points[i], points[j]));
        }
    }
    qsort(edges, numEdges, sizeof (Edge*), compareEdgesByLen);

    /* union-find: quickunion bilanciata */
    Node **partizione = malloc(sizeof (Node*) * n);
    for (int i=0; i<n; i++)
        partizione[i] = newNode(i);

    int edgesAdded = 0;
    int currEdge = 0;
    while (edgesAdded < n-1) {
        Edge *e = edges[currEdge++];
        Node *fatherA = NULL;
        Node *fatherB = NULL;
        /* find a */
        fatherA = partizione[e->a];
        while (fatherA->father != NULL) fatherA = fatherA->father;
        /* find b */
        fatherB = partizione[e->b];
        while (fatherB->father != NULL) fatherB = fatherB->father;
        
        if (fatherA->i != fatherB->i) {
            /* union */
            if (fatherB->height < fatherA->height)
                fatherB->father = fatherA;
            else if (fatherB->height > fatherA->height)
                fatherA->father = fatherB;
            else {
                fatherB->father = fatherA;
                fatherA->height += 1;
            }
            /* set matrix */
            adjacency[e->a * n + e->b] = true;
            adjacency[e->b * n + e->a] = true;
            edgesAdded++;
        }
    }

    for (int i=0; i<n; i++)
        free(partizione[i]);
    free(partizione);
    for (int i=0; i<numEdges; i++)
        free(edges[i]);
    free(edges);
}

bool areOverlapped(Point *p1, Point *p2) {
    return fabs(p1->x - p2->x) + fabs(p1->y - p2->y) < 0.0001;
}

double z(Point *x, Elem *elems) {
    double res = 0;
    for (Elem *i=elems; i!=NULL; i=i->next)
        res += i->w * dist(x, i->p);
    return res;
}

Point *coreLocate(Elem *anchors, double precision) {
    /* finding min anchor */
    Elem *bestAnchor = anchors;
    double bestVal = FLT_MAX;
    for (Elem *i=anchors; i!=NULL; i=i->next) {
        double valZ = z(i->p, anchors);
        if (valZ < bestVal) {
            bestVal = valZ;
            bestAnchor = i;
        }
    }

    double Rj[2] = {0,0}; /* vector not point */
    Point *q = bestAnchor->p;
    double wj = 0;
    for (Elem *i=anchors; i!=NULL; i=i->next) {
        double d = dist(q, i->p);
        if (d != 0) {
            Rj[0] += i->w / d * (q->x - i->p->x);
            Rj[1] += i->w / d * (q->y - i->p->y);
        } else wj += i->w;
    }
    double normRj = sqrt(pow(Rj[0], 2) + pow(Rj[1], 2));

    Point *xNew = newPoint(0,0);
    if (normRj <= wj) {
        xNew->x = q->x;
        xNew->y = q->y;
    }
    else { // x liberare la memoria in entrambi i casi
        Point *xOld = q;
        double L = 0;
        for (Elem *i=anchors; i!=NULL; i=i->next) {
            double d = dist(xOld, i->p);
            if (d != 0) L += i->w / d;
        }
        double fatt = (normRj - wj) / (L * normRj);
        xNew->x = xOld->x - (fatt * Rj[0]);
        xNew->y = xOld->y - (fatt * Rj[1]);

        Point *xTemp = newPoint(0,0);
        while (dist(xNew,xOld) > precision) {
            xOld = xNew;
            xNew = xTemp;
            xNew->x = 0; xNew->y = 0;

            double L = 0;
            for (Elem *i=anchors; i!=NULL; i=i->next) {
                Point *ai = i->p;
                double d = dist(xOld, ai);
                L += i->w / d;
                xNew->x += (i->w * ai->x) / d;
                xNew->y += (i->w * ai->y) / d;
            }
            xNew->x /= L;
            xNew->y /= L;
            xTemp = xOld;
        }
        free(xTemp);
    }
    do {
        Elem *i = anchors->next;
        free(anchors);
        anchors = i;
    } while (anchors != NULL);
    return xNew;
}
Point *locate(Point *terminals[], double weights[], int n, double precision) {
    /* n = # punti totali (terminals + sp) */
    /* removing anchors with weight 0*/
    Elem *anchors = NULL;
    for (int i=0; i<n; i++) {
        if (weights[i] != 0)
            anchors = newElem(terminals[i], weights[i], anchors);
    }
    return coreLocate(anchors, precision);
}
Point *simpleLocate(Point *terminals[], bool weights[], int n, double precision) {
    /* n = # punti totali (terminals + sp) */
    /* removing anchors with weight 0*/
    Elem *anchors = NULL;
    for (int i=0; i<n; i++) {
        if (weights[i] != 0)
            anchors = newElem(terminals[i], weights[i], anchors);
    }
    return coreLocate(anchors, precision);
}  // "simple" because weights are boolean

int closestPoint(int s, Point **p, int n) {
    /* restituisce l'indice del punto più vicino a quello di indice s */
    int closest = 0;
    double distance = FLT_MAX;
    for (int i=0; i<n; i++) {
        if (i != s) {
            double d = dist(p[s], p[i]);
            if (d < distance) {
                closest = i;
                distance = d;
            }
        }
    }
    return closest;
}


void printPoint(Point *p) {
    printf("(%.3f,%.3f)\n", p->x, p->y);
}
void printPointsPy(Point **p, int n) {
    for (int i=0; i<n; i++)
        printf("%.7f,%.7f\n", p[i]->x, p[i]->y);
}
void printPoints(Point **p, int n) {
    for (int i=0; i<n; i++)
        printf("(%.3f,%.3f)\n", p[i]->x, p[i]->y);
}
void printMatrixAdj(bool matrix[], int row, int columns) {
    for (int i=0; i<row; i++) {
        for (int j=0; j<columns; j++) {
            printf("%d ", matrix[columns*i +j]);
        }
        printf("\n");
    }
}
void printMatrixDouble(double matrix[], int row, int columns) {
    for (int i=0; i<row; i++) {
        for (int j=0; j<columns; j++) {
            printf("%.2f ", matrix[columns*i +j]);
        }
        printf("\n");
    }
}
void printMatrixAdjPy(bool matrix[], int row, int columns) {
    for (int i=0; i<row; i++) {
        for (int j=0; j<columns; j++) {
            printf("%d", matrix[columns*i +j]);
            if (j<columns-1) printf(" ");
        }
        printf("\n");
    }
}
void printMatrixDoublePy(double matrix[], int row, int columns) {
    for (int i=0; i<row; i++) {
        for (int j=0; j<columns; j++) {
            printf("%.0f", matrix[columns*i +j]);
            if (j<columns-1) printf(" ");
        }
        printf("\n");
    }
}
void printVariables(double x[], double y[], int n) {
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++)
            printf("%.3f ", x[s*n+t]);
        printf("  | ");
        for (int s2=0; s2<n-2; s2++)
            printf("%.3f ", y[s*(n-2)+s2]);
        printf("\n");
    }
}
void printVariablesPy(double x[], double y[], int n) {
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++) {
            printf("%.3f", x[s*n+t]);
            if (t<n-1) printf(" ");
        }
        printf("|");
        for (int s2=0; s2<n-2; s2++) {
            printf("%.0f", y[s*(n-2)+s2]);
            if (s2<n-3) printf(" ");
        }
        printf("\n");
    }
}
void sanityCheck(double x[], int n) {
    for (int t=0; t<n; t++) {
        bool found = false;
        for (int s=0; s<n-2; s++) {
            if (x[s*n+t] > 0.9) {
                if (found) printf("Terminal %d linked to more sp\n", t);
                found = true;
            }
        }
        if (!found) printf("Terminal %d not linked to any sp\n", t);
    }
}
int sanityCheckSilent(double x[], int n) {
    for (int t=0; t<n; t++) {
        bool found = false;
        for (int s=0; s<n-2; s++) {
            if (x[s*n+t] > 0.7) {
                if (found) return t;
                found = true;
            }
        }
        if (!found) return t;
    }
    return -1;
}
void printXSanityCheck(double x[], int n) {
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++) printf("%.3f ", x[s*n+t]);
        printf("\n");
    }
}
void printSolLP(double w[], int n) {
    // n = # terminali
    for (int s=0; s<n-2; s++) {
        int colIndex = 0;
        while (colIndex<n) {
            printf("%.3f", w[s*n*(n-2)+colIndex]);
            if (colIndex<n-1) printf(" ");
            colIndex++;
        }
        printf("|");
        for (int s2=0; s2<n-2; s2++) {
            printf("%.3f", w[s*n*(n-2)+colIndex]);
            if (s2<n-3) printf(" ");
            colIndex++;
        }
        printf("\n");
    }
}
void printSolLP_v2(double w[], int n) {
    // n = # terminali
    int numCol = 0;
    int numRow = 0;
    int index = 0;
    for (int i=0; i<(n-2)*(2*n-2); i++) {
        if (numCol - n == numRow) printf("%.3f ", 0.0);
        else printf("%.3f ", w[index++]);
        numCol++;
        if (numCol == n) printf(" | ");
        if (numCol == 2*n-2) {
            numCol = 0;
            numRow++;
            printf("\n");
        }
    }
}
void printSolLP_v3(double x[], double w[], int n) {
    // n = # terminali
    int index = 0;
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++) {
            printf("%.3f ", x[s*n + t]);
        }
        printf("| ");
        for (int r=0; r<n-2; r++) {
            double elem = (r != s) ? w[index++] : 0;
            printf("%.3f ", elem);
        }
        printf("\n");
    }
}

#define INITIAL_TEMPERATURE     10
#define SF                       0.9  
#define MAX_NUM_RIPETIZIONI     40
#define FINAL_SATURATION         0.999
#define PRECISION_LOCATE_MID     0.001
#define PRECISION_LOCATE_FINAL   0.0000000001
#define SPOSTAMENTO_MEDIO_MIN    0.005
#define SPOSTAMENTO_MEDIO_MINF   0.0001

#define DEG_MAX_SURPLUS          0.1
#define CLOSE_DISTANCE           0.01

typedef enum {
    ALL,
    ROW
} normalization;

int main(int argc, char **argv) {
    bool debug = false;
    bool verbose = false;
    normalization actualNormalization = ROW;
    if (argc == 3 && strcmp("-n", argv[1]) == 0) {
        if (strcmp("all", argv[2]) == 0) actualNormalization = ALL;
        else if (strcmp("row", argv[2]) == 0) actualNormalization = ROW;
    }

    int n = 0;
    scanf(" %d\n", &n);
    Point **a = malloc(sizeof (Point*) * (n)); /* anchors/terminals */

    double xmin = FLT_MAX;
    double ymin = FLT_MAX;
    double xmax = -FLT_MAX;
    double ymax = -FLT_MAX;
    for (int i=0; i<n; i++) {
        a[i] = newPoint(0,0);
        scanf(" %lf %lf\n", &(a[i]->x),&(a[i]->y));
        /* determinazione len diagonale rettangolo */
        if (a[i]->x < xmin) xmin = a[i]->x;
        if (a[i]->y < ymin) ymin = a[i]->y;
        if (a[i]->x > xmax) xmax = a[i]->x;
        if (a[i]->y > ymax) ymax = a[i]->y;
    }
    double T = sqrt(pow(xmax - xmin, 2) + pow(ymax - ymin, 2)) * 10;
    fprintf(stderr, "T = %f\n", T);
    Point **sp = malloc(sizeof (Point*) * (n-2));

    double *x = malloc(sizeof (double) * n * (n-2));
    double *y = malloc(sizeof (double) * (n-2) * (n-2));
    
    long int seed = time(NULL);
    srand(seed);
    if (verbose) fprintf(stderr, "Seed = %ld\n", seed);
    /* disponi sp a caso */
    for (int s=0;s<n-2;s++)
        sp[s] = newPoint((double)rand()/RAND_MAX * xmax, (double)rand()/RAND_MAX * ymax);

    const double delta = 700; /* max esponente e^ */
    double *esponenti = malloc(sizeof (double) * (n-2)); /* per il calcolo di x */
    Point **newSp = malloc(sizeof (Point*) * (n-2)); /* per locazione nuovi sp */
    double satX = 0;
    double satY = 0;
    double spostamentoMedio = 0;
    double distanzaMedia = 0;

    int numFrames = 0; /* x animazione */

    bool start = false;
    bool alreadyRaised = false;

    /* LP problem to solve with GLPK */
    double *w = malloc(sizeof (double) * (n-2) * (n - 3));
    int ne = (n-2)*(n-3) + 4*(n-2)*(n-3) + (n-2)*(n-3) + (n-2)*(n-3)*((int) pow(2,n-4));
    int *ia = malloc(sizeof (int) * (ne + 1));
    int *ja = malloc(sizeof (int) * (ne + 1));
    double *ar = malloc(sizeof (double) * (ne + 1));
    glp_prob *lp;
    Point **allAnchors = malloc(sizeof (Point*) * (2*n-2)); // per weiszfeld
    double *pesiFinali = malloc(sizeof (double) * (n-2) * (2*n-2));

    double *esponentiX = malloc(sizeof (double) * (n-2)); // per nuovo calcolo x
    double *allEsponentiY = malloc(sizeof (double) * (n-2) * (n-3)); // per nuovo calcolo y
    // n-3 elementi per riga xk le y sulla diagonale le metto a 0 direttamente


    double satW = 0;
    double saturation = 0;
    bool *movable = malloc(sizeof (bool) * (n-2));
    /* start_annelleaing: */
    while (saturation < 0.99) {
        int counter = 0;
        do {
            if (debug) printPointsPy(sp, n-2);
            int numAlreadyRepeated = 0;
            int wrongDegS = 0;

            /* aggiorno movable */
            for (int s=0; s<n-2; s++) {
                int closest = closestPoint(s, sp, n-2);
                double degS = 0;
                for (int t=0; t<n; t++) degS += x[s*n + t];
                double degC = 0;
                for (int t=0; t<n; t++) degC += x[closest*n + t];
                if (dist(sp[s], sp[closest]) < CLOSE_DISTANCE && degS + degC > 2) {
                    movable[s] = false;
                    movable[closest] = false;
                } else movable[s] = true;
            }
            if (debug) { /* print movable */
                printf("[ ");
                for (int s=0; s<n-2; s++) printf("%d ", movable[s]);
                printf("]\n");
            }
            
            while (wrongDegS != -1) {
                /* calcolo x */
                satX = 0;
                for (int t=0; t<n; t++) {
                    double minEsp = FLT_MAX;
                    double sum = 0;
                    for (int s=0; s<n-2; s++) {
                        esponentiX[s] = dist(a[t], sp[s]) / T;
                        if (esponentiX[s] < minEsp) minEsp = esponentiX[s];
                        x[s*n+t] = exp(-esponentiX[s]);
                        sum += x[s*n+t];
                    }
                    if (sum == 0) {
                        for (int s=0; s<n-2; s++) {
                            double sigma = dist(sp[s], a[t]) / T;
                            if (sigma <= minEsp + delta) {
                                double sum = 0;
                                for (int i=0; i<(n-2); i++) sum += exp(sigma-esponentiX[i]);
                                x[s*n + t] = (1.0 / sum);
                            } else x[s*n + t] = 0;
                            double xij = x[s*n+t];
                            satX += (xij < 1.0/(n-2)) ? 1 - (n-2)*xij : (1.0/(n-3))*((n-2)*xij-1);
                        }
                    } else {
                        for (int s=0; s<n-2; s++) {
                            x[s*n+t] /= sum;
                            double xij = x[s*n+t];
                            satX += (xij < 1.0/(n-2)) ? 1 - (n-2)*xij : (1.0/(n-3))*((n-2)*xij-1);
                        }
                    }
                }

                /* rilocalizzazione steiner points: */
                /* se sp collegato ai terminali per grado tot > 2 sposto il sp */
                /* meno collegato ai terminali lì (coincidono) */
                /*(else il problema di LP è insoddisfacibile (vincolo grado3 e connettività))*/
                /*una volta che ho spostato un punto sopra un altro non li posso più muovere,*/
                /*else quel punto ritorna ad avere grado sbagliato*/
                /*(per questo mantengo array movable)*/

                /* riazzero movable al più 1 volte */
                /* else potrebbe non riuscire a portare tutti a grado < 2 */
                double minDeg = FLT_MAX;
                double maxDeg = -FLT_MAX;
                int minDegS = -1; // indice del sp *spostabile* con grado minore
                wrongDegS = -1;
                bool smallExcess = false;
                for (int s=0; s<n-2; s++) {
                    double deg = 0;
                    for (int t=0; t<n; t++) deg += x[s*n + t];
                    if (debug) printf("%.2f ", deg);
                    if (deg < minDeg && movable[s]) {
                        minDeg = deg;
                        minDegS = s;
                    }
                    if (deg > 2 && deg > maxDeg) {
                        maxDeg = deg;
                        wrongDegS = s;
                        if (deg-2 < DEG_MAX_SURPLUS) { // keep trying to add noise
                            if (debug) printf("surplus %d = %g\n", s, deg-2);
                            sp[s]->x = sp[s]->x + (((double) rand()/RAND_MAX) / 20 - 0.025);
                            sp[s]->y = sp[s]->y + (((double) rand()/RAND_MAX) / 20 - 0.025);
                            smallExcess = true;
                        }
                    }
                }
                if (debug) printf("\n");
                if (smallExcess) continue;
                if (wrongDegS != -1) {
                    if (minDegS == -1 && numAlreadyRepeated < 1) {
                        /* per il secondo tentativo rendo tutti movable */
                        for (int s=0; s<n-2; s++) movable[s] = true;
                        numAlreadyRepeated++;
                        wrongDegS = 0;
                        fprintf(stderr, "\n\nsecond try\n----------\n");
                        continue;
                    }
                    if (minDegS == -1 && numAlreadyRepeated == 1) {printf("INFINITE LOOP wrongDeg  [seed = %ld]\n", seed); return -1;}
                    if (debug) printf("Rilocated: moved %d (%.3f,%.3f) to %d (%.3f,%.3f)\n", minDegS, sp[minDegS]->x, sp[minDegS]->y, wrongDegS, sp[wrongDegS]->x, sp[wrongDegS]->y);
                    sp[minDegS]->x = sp[wrongDegS]->x;
                    sp[minDegS]->y = sp[wrongDegS]->y;
                    movable[minDegS] = false;
                    movable[wrongDegS] = false;
                } else {
                    for (int s=0; s<n-2; s++) movable[s] = true; // reset
                    /* calcolo y */
                    if (actualNormalization == ALL) {
                        satY = 0;
                        double minEsp = FLT_MAX;
                        double sum = 0;
                        for (int s1=0; s1<n-2; s1++) {
                            int currIndex = 0;
                            for (int s2=0; s2<n-2; s2++) {
                                if (s2 != s1) {
                                    allEsponentiY[s1*(n-3) + currIndex] = dist(sp[s1], sp[s2]) / (2*T);
                                    if (allEsponentiY[s1*(n-3) + currIndex] < minEsp) minEsp = allEsponentiY[s1*(n-3) + currIndex];
                                    y[s1*(n-2) + s2] = exp(-allEsponentiY[s1*(n-3) + currIndex]);
                                    sum += y[s1*(n-2)+s2];
                                    currIndex++;
                                }
                                else y[s1*(n-2)+s2] = 0;
                            }
                        }
                        if (sum == 0) {
                            for (int s1=0; s1<n-2; s1++) {
                                for (int s2=0; s2<n-2; s2++) {
                                    if (s2 != s1) {
                                        double sigma = dist(sp[s1], sp[s2]) / (2*T);
                                        if (sigma <= minEsp + delta) {
                                            double sum = 0;
                                            for (int i=0; i<(n-2)*(n-3); i++) sum += exp(sigma-allEsponentiY[i]);
                                            y[s1*(n-2) + s2] = (1.0 / sum) * (2*n-6);
                                        } else y[s1*(n-2) + s2] = 0;
                                    }
                                    double yij = y[s1*(n-2)+s2];
                                    satY += (yij < 1.0/((n-2)*(n-2))) ? 1 - ((n-2)*(n-2))*yij : (1.0/((n-2)*(n-2)-1))*(((n-2)*(n-2))*yij-1);
                                }
                            }
                        } else {
                            for (int s1=0; s1<n-2; s1++) {
                                for (int s2=0; s2<n-2; s2++) {
                                    y[s1*(n-2)+s2] =  (y[s1*(n-2)+s2] / sum) * (2*n-6); // collegamenti steiner-steiner contati 2 volte
                                    double yij = y[s1*(n-2)+s2];
                                    satY += (yij < 1.0/((n-2)*(n-2))) ? 1 - ((n-2)*(n-2))*yij : (1.0/((n-2)*(n-2)-1))*(((n-2)*(n-2))*yij-1);
                                }
                            }
                        }
    

                    } else if (actualNormalization == ROW) {
                        satY = 0;
                        double minEsp = FLT_MAX;
                        for (int s1=0; s1<n-2; s1++) {
                            double sum = 0;
                            int currIndex = 0;
                            double sumX = 0;
                            for (int t=0; t<n; t++) sumX += x[s1*n + t];
                            for (int s2=0; s2<n-2; s2++) {
                                if (s2 != s1) {
                                    allEsponentiY[s1*(n-3) + currIndex] = dist(sp[s1], sp[s2]) / (2*T);
                                    if (allEsponentiY[s1*(n-3) + currIndex] < minEsp) minEsp = allEsponentiY[s1*(n-3) + currIndex];
                                    y[s1*(n-2) + s2] = exp(-allEsponentiY[s1*(n-3) + currIndex]);
                                    sum += y[s1*(n-2)+s2];
                                    currIndex++;
                                }
                                else y[s1*(n-2)+s2] = 0;
                            }
                            if (sum == 0) {
                                for (int s2=0; s2<n-2; s2++) {
                                    if (s2 != s1) {
                                        double sigma = dist(sp[s1], sp[s2]) / (2*T);
                                        if (sigma <= minEsp + delta) {
                                            double sum = 0;
                                            for (int i=0; i<(n-3); i++) sum += exp(sigma-allEsponentiY[s1*(n-3) + i]);
                                            y[s1*(n-2) + s2] = (1.0 / sum) * (3 - sumX);
                                        } else y[s1*(n-2) + s2] = 0;
                                    }
                                    double yij = y[s1*(n-2)+s2];
                                    satY += (yij < 1.0/(n-2)) ? 1 - (n-2)*yij : (1.0/(n-2-1))*((n-2)*yij-1);
                                }
                            } else {
                                for (int s2=0; s2<n-2; s2++) {
                                    y[s1*(n-2)+s2] =  (y[s1*(n-2)+s2] / sum) * (3 - sumX); // collegamenti steiner-steiner contati 2 volte
                                    double yij = y[s1*(n-2)+s2];
                                    satY += (yij < 1.0/(n-2)) ? 1 - (n-2)*yij : (1.0/(n-2-1))*((n-2)*yij-1);
                                }
                            }
                        }            
                    }
                }
            }

            satX /= (n-2) * (n);
            satY /= (n-2) * (n-2);

            /* calcolo dei pesi: GLPK */
            /*  Y ha n-3 colonne */
            lp = glp_create_prob();
            glp_set_obj_dir(lp, GLP_MIN);

            glp_add_rows(lp, n-2 + 2*(n-2)*(n-3) + ((n-2)*(n-3))/2 + ((int) pow(2,n-2) -2));
            int rowsAdded = 0;
            
            for (int i=0; i<n-2; i++) { // grado 3 per i punti di Steiner
                double sumX = 0;
                for (int t=0; t<n; t++) sumX += x[i*n + t];
                if (sumX > 3) printf("ERROR: sumX > 3 (%.3f)\n", sumX);
                glp_set_row_bnds(lp, ++rowsAdded, GLP_FX, 3.0 - sumX, 3.0 - sumX);
            }

            glp_add_cols(lp, (n-2)*(n-3)+1);
            for (int i=0; i<(n-2)*(n-3); i++) { // w_ij (senza diagonale Y)
                glp_set_col_bnds(lp, i+1 , GLP_DB, 0.0, 1.0);
                glp_set_obj_coef(lp, i+1, 0.0);
            }
            // delta (max)
            glp_set_col_bnds(lp, (n-2)*(n-3)+1 , GLP_FR, 0.0, 0.0);
            glp_set_obj_coef(lp, (n-2)*(n-3)+1, 1.0);


            /*** imposta matrice ***/
            int numCol = n-3;
            /*   dispongo i w_ij concatenando le righe */
            //   - grado 1 per i terminali
            int nextElemIndex = 1;

            /*   1   2       f   f+1 f+2 ...         (f = numCol) */
            /* w11 w12 ... w1f   w21 w22 ... w2f ... */

            //   - grado 3 per i punti di Steiner
            for (int i=0; i<n-2; i++) {
                for (int j=1; j<=n-3; j++) {
                    ia[nextElemIndex] = i+1;
                    ja[nextElemIndex] = i*numCol + j;
                    ar[nextElemIndex] = 1.0;
                    nextElemIndex++;
                }
            }

            int numNextRiga = n-2 + 1;

            //   - delta (max abs)
            int shift = 0; // xk da Y ho rimosso diagonale                                                 
            for (int i=0; i<(n-2)*(n-3); i++) {
                ia[nextElemIndex] = numNextRiga;
                ja[nextElemIndex] = i+1;
                ar[nextElemIndex] = -1.0;
                nextElemIndex++;
            
                int originalRow = i / numCol;                                                              
                int originalCol = i % numCol;                                                              
                if (originalCol == 0) shift = 0;
                double y_ij = 0;
                if (originalCol == originalRow) shift = 1;
                y_ij = y[(n-2)*originalRow + originalCol + shift];                                         
                
                
                ia[nextElemIndex] = numNextRiga;
                ja[nextElemIndex] = (n-2)*(n-3)+1;
                ar[nextElemIndex] = 1.0;
                glp_set_row_bnds(lp, ++rowsAdded, GLP_LO, -y_ij, 0.0);
                numNextRiga++;
                nextElemIndex++;
                
                
                ia[nextElemIndex] = numNextRiga;
                ja[nextElemIndex] = i+1;
                ar[nextElemIndex] = 1.0;
                nextElemIndex++;
            
                ia[nextElemIndex] = numNextRiga;
                ja[nextElemIndex] = (n-2)*(n-3)+1;
                ar[nextElemIndex] = 1.0;
                glp_set_row_bnds(lp, ++rowsAdded, GLP_LO, y_ij, 0.0);
                numNextRiga++;
                nextElemIndex++;
            }
            /**/



            for (int i=0; i<((n-2)*(n-3))/2; i++) { // simmetria Y
                glp_set_row_bnds(lp, ++rowsAdded, GLP_FX, 0.0, 0.0);
            }
            //   - simmetria Y
            for (int i=1; i<=n-3; i++) {
                for (int j=i; j<=n-3; j++) {
                    // printf("%d,%d  %d -- %d\n", i,j,(i-1)*numCol + n + j,j*numCol + n + i);
                    ia[nextElemIndex] = numNextRiga;
                    ja[nextElemIndex] = (i-1)*numCol + j;
                    ar[nextElemIndex] = 1.0;
                    nextElemIndex++;

                    ia[nextElemIndex] = numNextRiga;
                    ja[nextElemIndex] = j*numCol + i;
                    ar[nextElemIndex] = -1.0;
                    nextElemIndex++;
                    numNextRiga++;
                }
            }

            for (int i=0; i<(int) pow(2,n-2)-2; i++) { // connettività
                glp_set_row_bnds(lp, ++rowsAdded, GLP_LO, 1.0, 1.0);
            }
            //   - connettività
            for (int i=1; i<=(int) pow(2,n-2)-2; i++) {
                for (int s1=0; s1<n-2; s1++) {
                    if ((int) (i/pow(2,s1)) % 2 == 1) {
                        for (int s2=0; s2<n-2; s2++) {
                            if ((int) (i/pow(2,s2)) % 2 == 0) {
                                ia[nextElemIndex] = numNextRiga;
                                ja[nextElemIndex] = s1*numCol + (s2<=s1 ? s2+1 : s2);
                                ar[nextElemIndex] = 1.0;
                                nextElemIndex++;
                            }
                        }
                    }
                }
                numNextRiga++;
            }

            glp_load_matrix(lp, ne, ia, ja, ar);
            /* per risolvere con dual simplex *
             glp_smcp parm;
             glp_init_smcp(&parm);
             parm.meth = GLP_DUAL;
             glp_simplex(lp, &parm);
            /**/
            /* per disattivare output */
             glp_smcp parm;
             glp_init_smcp(&parm);
             parm.msg_lev = GLP_MSG_OFF;
             glp_simplex(lp, &parm);
            /**/
            // glp_simplex(lp, NULL);

            satW = 0;
            for(int i=0; i<(n-2)*(n-3); i++) {
                w[i] = glp_get_col_prim(lp, i+1);
                satW += (w[i] < 2.0/(n-2)) ? 1 - (w[i]*(n-2) / 2) : ((w[i]*(n-2) - 2) / 2);
            }

            glp_delete_prob(lp);


            satW /= (n-2) * (n-3);
            saturation = (satX + satW) / 2;

            /* locate sp */
            int indexW = 0;
            for (int s=0; s<n-2; s++) {
                for (int t=0; t<n; t++) {
                    pesiFinali[s*(2*n-2) + t] = x[s*n + t];
                }
                for (int r=0; r<n-2; r++) {
                    if (s != n-3 || r != n-3)
                        pesiFinali[s*(2*n-2) + n + r] = (s == r) ? 0 : w[indexW++];
                    else pesiFinali[s*(2*n-2) + n + r] = 0;
                    // else problema con ultimo elemento diagonale (non posso prendere indexW)
                }
            }

            spostamentoMedio = 0;
            for (int t=0; t<n; t++) allAnchors[t] = a[t];
            for (int s=0; s<n-2; s++) allAnchors[n+s] = sp[s];
            for (int s=0; s<n-2; s++) {
                newSp[s] = locate(allAnchors, &pesiFinali[s*(2*n-2)], 2*n-2, PRECISION_LOCATE_MID);
                spostamentoMedio += dist(newSp[s], sp[s]);
                // aggiungo il rumore solo ai punti che si sovrapporrebbero
                if (saturation < 0.99) {
                    for (int s2=0; s2<s; s2++) {
                        if (areOverlapped(newSp[s], newSp[s2])) {
                            newSp[s]->x += ((double) rand()/RAND_MAX) / 20 - 0.025;
                            newSp[s]->y += ((double) rand()/RAND_MAX) / 20 - 0.025;

                            movable[s] = false;
                            movable[s2] = false;
                            break;
                        }
                    }
                }
            }
            spostamentoMedio /= (n-2);
            for (int s=0; s<n-2; s++) {
                free(sp[s]);
                sp[s] = newSp[s];
            }

            distanzaMedia = 0;
            for (int s1=0; s1<n-2;s1++)
                for (int s2=s1+1; s2<n-2; s2++) distanzaMedia += dist(sp[s1], sp[s2]);
            distanzaMedia /= (n-2)*(n-3)/2;

            /**/
            if (verbose && T < 0.3) {
                fprintf(stdout, "Satx = %.3f\n", satX);
                fprintf(stdout, "Satw = %.3f\n", satW);
                fprintf(stdout, "saturation = %.3f\n", saturation);
                printVariables(x,y,n);
                printf("w:\n");
                // for (int i=0; i<(n-2)*(2*n-2); i++) printf("%.2f\n", w[i]);
                printSolLP_v3(x, w, n);
                    // printf("Spostamento medio = %.3f\n", spostamentoMedio);
                    // printf("Distanza media    = %.3f\n", distanzaMedia);
                    // printf("Steiners:\n");
                printPoints(sp,n-2);
                // if (counter+1 == ripetizioniMax) fprintf(stderr, "end cause max repetition\n");
                printf("---------------------------------------------------\n");
            }
            /**/


            //satX = 1000;
        } while ( ++counter < MAX_NUM_RIPETIZIONI && spostamentoMedio > SPOSTAMENTO_MEDIO_MIN ); /* pos no cambiata molto */
        T *= SF;
        numFrames++;

        /**/
        if (verbose) {
            printVariables(x,y,n);
            // printf("\n");
            printPoints(sp,n-2);
            fprintf(stdout, "Satx = %.3f\n", satX);
            fprintf(stdout, "\n\n================\n T = %.7f\n================\n", T);
            // printf("Satx = %.3f\n", satX);
            // printf("\n\n================\n T = %.7f\n================\n", T);
        }
        /**/
        // fprintf(stderr, "\n\n================\n T = %.7f\n================\n", T);

        /* stampe per animazione *
        // printf(" %2d\n====\n", numFrames);
        printVariablesPy(x,y,n);
        printf("\n");
        printPointsPy(sp,n-2);
        printf("\n");
        /**/

    }

    /* rimuovi punti di steiner sovrapposti */
    int numSp = n-2;
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++) {
            if (areOverlapped(a[t], sp[s])) {
                sp[s]->x = NAN; // cancellazione logica
                numSp--;
                break;
            }
        }
    }

    Point **points = malloc(sizeof (Point*) * (n + numSp));
    for (int t=0; t<n; t++) points[t] = a[t];
    int nextIndex = n;
    for (int s=0; s<n-2; s++) if (!isnan(sp[s]->x)) points[nextIndex++] = sp[s];

    bool *adjacency = malloc(sizeof (bool) * (n + numSp) * (n + numSp));
    emst(points, n+numSp, adjacency);

    /* cancella sp con grado 2 (possono comparire dopo che calcolato emst) */
    int removed = 0;
    bool removedOne = true;
    while (removedOne) {
        removedOne = false;
        for (int s=0; s<numSp; s++) {
            int degree = 0;
            int j1 = -1;
            int j2 = -1;
            for (int j=0; j<n+numSp; j++) {
                if (adjacency[(n+s)*(n+numSp) + j] == 1) {
                    degree++;
                    if (j1 == -1) j1 = j;
                    else j2 = j;
                }
            }
            if (degree < 3) {
                removedOne = true;
                free(points[n+s]);
                points[n+s] = points[n+numSp-1];
                numSp--;
                s--;
                emst(points, n+numSp, adjacency);
            }
        }
    }

    /* steinerization */
    Point **nuoviSp = malloc(sizeof (Point*) * n); // non dovrei aggiungerne più di n
    int spAdded = 0;
    SteinerItem **updated = malloc(sizeof (SteinerItem*) * n);
    for (int t=0; t<n; t++) updated[t] = NULL;
    for (int t=0; t<n; t++) { // controllo solo i terminali
        // con i sp c'era il problema che a volte angolo vicino a 120° ma minore
        for (int j1=0; j1<n+numSp; j1++) {
            if (adjacency[(n+numSp)*t+j1]==1) {
                SteinerItem *si1 = updated[t];
                while (si1 != NULL && si1->index != j1) si1 = si1->next;
                Point *point1 = si1 != NULL ? si1->sp : points[j1];
                for (int j2=j1+1; j2<n+numSp; j2++) {
                    if (adjacency[(n+numSp)*t+j2]==1) {
                        if (verbose) fprintf(stderr, "%2d -- %2d -- %2d\n", t,j1,j2);
                        SteinerItem *si2 = updated[t];
                        while (si2 != NULL && si2->index != j2) si2 = si2->next;
                        Point *point2 = si2 != NULL ? si2->sp : points[j2];
                        double ax = point1->x - points[t]->x;
                        double ay = point1->y - points[t]->y;
                        double bx = point2->x - points[t]->x;
                        double by = point2->y - points[t]->y;
                        double cosAlpha = (ax*bx + ay*by) / (dist(points[t],point1)*dist(points[t],point2));
                        if (cosAlpha > -0.49) { // cos(120°) = -0.5
                            Point *anchors[3] = {points[t], point1, point2};
                            bool weights[3] = {1,1,1};
                            if (spAdded == n) fprintf(stderr, "Added too many sps!!!\n");
                            else {
                                nuoviSp[spAdded++] = simpleLocate(anchors, weights, 3, PRECISION_LOCATE_MID);
                                // segno cambiamento per gli altri due punti se sono terminali
                                if (si1 == NULL && j1<n)
                                    updated[j1] = newItem(t, nuoviSp[spAdded-1], updated[j1]);
                                if (si2 == NULL && j2<n)
                                    updated[j2] = newItem(t, nuoviSp[spAdded-1], updated[j2]);
                                if (verbose) fprintf(stderr, "%2d -- %2d -- %2d  --> ", t,j1,j2);
                                if (verbose) fprintf(stderr, "(%.4f,%.4f)\n", nuoviSp[spAdded-1]->x,nuoviSp[spAdded-1]->y);
                            }
                        }
                    }
                }
            }
        }

    }
    // free tutti i SteinerItem
    for (int t=0; t<n; t++) {
        SteinerItem *si = updated[t];
        while (si != NULL) {
            SteinerItem *curr = si;
            si = si->next;
            free(curr);
        }
    }
    free(updated);

    
    // - aggiorno i points x ricalcolare emst
    Point **finalPpoints = malloc(sizeof (Point*) * (n + numSp + spAdded));
    for (int i=0; i<n+numSp; i++) finalPpoints[i] = points[i];
    for (int i=0; i<spAdded; i++) {
        finalPpoints[n+numSp+i] = nuoviSp[i];
    }
    free(adjacency);
    adjacency = malloc(sizeof (bool) * (n + numSp + spAdded) * (n + numSp + spAdded));
    numSp += spAdded;
    free(points);
    emst(finalPpoints, n+numSp, adjacency);

    Point **tmpPoints = malloc(sizeof (Point*) * numSp);
    spostamentoMedio = 1;
    int numIteration = 0;
    while (spostamentoMedio > SPOSTAMENTO_MEDIO_MINF && numIteration < MAX_NUM_RIPETIZIONI) {
        spostamentoMedio = 0;
        numIteration++;
        for (int s=0; s<numSp; s++) {
            tmpPoints[s] = simpleLocate(finalPpoints, &adjacency[(n+numSp)*(n+s)], n+numSp, PRECISION_LOCATE_FINAL);
            spostamentoMedio += dist(tmpPoints[s], finalPpoints[n+s]);
        }
        spostamentoMedio /= (numSp);
        for (int s=0; s<numSp; s++) {
            free(finalPpoints[n+s]);
            finalPpoints[n+s] = tmpPoints[s];
        }
    }
    /**/

    /**/

    /* restituisco i sp */
    for (int s=0; s<numSp; s++) {
        printf("%.10f,%.10f\n", finalPpoints[n+s]->x, finalPpoints[n+s]->y);
    }

    return 0;
}