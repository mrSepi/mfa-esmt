#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
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
Point *locate(Point *terminals[], double weightsT[], Point *sp[], bool weightsS[], int n, double precision) {
    /* n = # punti totali (terminals + sp) */
    /* removing anchors with weight 0*/
    int numT = (n+2)/2;
    int numS = n-numT;
    Elem *anchors = NULL;
    for (int i=0; i<numT; i++) {
        if (weightsT[i] != 0)
            anchors = newElem(terminals[i], weightsT[i], anchors);
        else n--; /* true number of anchors */
    }
    for (int i=0; i<numS; i++) {
        if (weightsS[i] != 0)
            anchors = newElem(sp[i], weightsS[i], anchors);
        else n--; /* true number of anchors */
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
            printf("%.0f ", matrix[columns*i +j]);
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
void printVariables(double x[], bool y[], int n) {
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++)
            printf("%.2f ", x[s*n+t]);
        printf("  | ");
        for (int s2=0; s2<n-2; s2++)
            printf("%d ", y[s*(n-2)+s2]);
        printf("\n");
    }
}
void printVariablesPy(double x[], bool y[], int n) {
    for (int s=0; s<n-2; s++) {
        for (int t=0; t<n; t++) {
            printf("%.3f", x[s*n+t]);
            if (t<n-1) printf(" ");
        }
        printf("|");
        for (int s2=0; s2<n-2; s2++) {
            printf("%d", y[s*(n-2)+s2]);
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


#define INITIAL_TEMPERATURE     10
#define SF                       0.9  
#define MAX_NUM_RIPETIZIONI     40
#define FINAL_SATURATION         0.999
#define PRECISION_LOCATE_MID     0.001
#define PRECISION_LOCATE_FINAL   0.0000000001
#define SPOSTAMENTO_MEDIO_MIN    0.005
#define SPOSTAMENTO_MEDIO_MINF   0.0001

#define SAT_THRESHOLD            0.8
#define UB_DEG                   3.6
#define LB_DEG                   2.5

int main() {
    bool verbose = false;

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
    Point **sp = malloc(sizeof (Point*) * (n-2));

    double *x = malloc(sizeof (double) * n * (n-2));
    bool   *y = malloc(sizeof (bool) * (n-2) * (n-2));
    

    long int seed = time(NULL);
    srand(seed);
    if (verbose) fprintf(stderr, "Seed = %ld\n", seed);
    /* disponi sp a caso */
    for (int s=0;s<n-2;s++)
        sp[s] = newPoint((double)rand()/RAND_MAX * xmax, (double)rand()/RAND_MAX * ymax);

    const double delta = 700; /* max esponente e^ */ /* TODO no magic numbers */
    double *esponenti = malloc(sizeof (double) * (n-2)); /* per il calcolo di x */
    Point **newSp = malloc(sizeof (Point*) * (n-2)); /* per locazione nuovi sp */
    double satX = 0;
    double spostamentoMedio = 0;
    double distanzaMedia = 0;

    int numFrames = 0; /* x animazione */

    bool start = false;
    double *degree = malloc(sizeof (double) * (n-2)); /* il grado di ogni sp (deve essere 3)*/
    bool alreadyRaised = false;

    /* start_annelleaing */
    while (satX < FINAL_SATURATION) {
        int counter = 0;
        do {
            /* prima di calcolare x, quando satX inizia ad avvicinarsi a 1, sposto i sp */
            /* ~~in modo che ogni sp abbia grado 3~~ */
            /* dalla prima volta che satX > 0.8 */
            if (!start && satX>SAT_THRESHOLD) start = true;
            if (start) {
                for (int s=0; s<n-2; s++) {
                    degree[s] = 0;
                    for (int t=0; t<n; t++) degree[s] += x[s*n+t];
                    for (int s2=0; s2<n-2; s2++) degree[s] += y[s*(n-2)+s2];
                }
                for (int s=0; s<n-2; s++) {
                    if (degree[s] < LB_DEG) { // check threshold, no faccio distinzione tra deg 2 e 1
                        int s2 = 0;
                        while (s2 < n-2 && degree[s2] <= UB_DEG) s2++;
                        if (s2 < n-2) {
                            /* sposto s vicino a s2 */
                            sp[s]->x = sp[s2]->x;
                            sp[s]->y = sp[s2]->y;
                            sp[s]->x  += ((double) rand()/RAND_MAX) / 20 - 0.025;
                            sp[s]->y  += ((double) rand()/RAND_MAX) / 20 - 0.025;               
                            degree[s2]--;
                        }
                    }
                }
            }

            /* calcolo x */
            satX = 0;
            for (int t=0; t<n; t++) {
                double sum = 0;
                for (int s=0; s<n-2; s++) {
                    x[s*n+t] = exp(-dist(a[t], sp[s]) / T);
                    sum += x[s*n+t];
                }
                if (sum == 0) {
                    double maxEsp = -FLT_MAX;
                    for (int s=0;s<n-2;s++) {
                        esponenti[s] = -dist(sp[s], a[t]) / T;
                        if (esponenti[s] > maxEsp) maxEsp = esponenti[s];
                    }
                    for (int s=0; s<n-2; s++) {
                        if (esponenti[s] > maxEsp - delta) {
                            double sum = 0;
                            for (int s2=0;s2<n-2;s2++) {
                                sum += exp(esponenti[s2]-esponenti[s]);
                            }
                            x[s*n+t] = 1.0/sum;
                        } else x[s*n+t] = 0;
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
            satX /= n * (n-2);

            /* calcolo y (emst) */
            emst(sp, n-2, y);

            /* locate sp */
            spostamentoMedio = 0;
            for (int s=0; s<n-2; s++) {
                newSp[s] = locate(a, &x[s*n], sp, &y[s*(n-2)], 2*n-2, PRECISION_LOCATE_MID);
                spostamentoMedio += dist(newSp[s], sp[s]);
                // aggiungo il rumore solo ai punti che si sovrapporrebbero
                if (satX < FINAL_SATURATION) {
                    for (int s2=0; s2<s; s2++) {
                        if (areOverlapped(newSp[s], newSp[s2])) {
                            newSp[s]->x += ((double) rand()/RAND_MAX) / 20 - 0.025;
                            newSp[s]->y += ((double) rand()/RAND_MAX) / 20 - 0.025;
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

        } while ( ++counter < MAX_NUM_RIPETIZIONI && spostamentoMedio > SPOSTAMENTO_MEDIO_MIN ); /* pos no cambiata molto */
        T *= SF;
        numFrames++;

        /* stampe per animazione *
        // printf(" %2d\n====\n", numFrames);
        printVariablesPy(x,y,n);
        printf("\n");
        printPointsPy(sp,n-2);
        printf("\n");
        /**/
    }

    //printVariablesPy(x,y,n);
    //printf("\n");
    //printPointsPy(sp,n-2);

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

    /* restituisco i sp */
    for (int s=0; s<numSp; s++) {
        printf("%.10f,%.10f\n", finalPpoints[n+s]->x, finalPpoints[n+s]->y);
    }

    return 0;
}