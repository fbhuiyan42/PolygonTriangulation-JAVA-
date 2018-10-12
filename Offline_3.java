/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package offline_3;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;


class Edge implements Comparable<Edge>
{
    Point p;
    Point q;

    public int compareTo(Edge t)
    {
       if(t.p.y == t.q.y)
            {
                if(p.y == q.y)
                {
                    if(p.y < t.p.y) return 1;
                    else return 0;
                }
                if(Point.convex(t.p,q,p)) return 1;
                else return 0;
            }
            else if(p.y == q.y)
            {
                if(Point.convex(p,t.q,t.p)) return 0;
                else return 1;
            }
            else if(p.x < t.p.x)
            {
                if(Point.convex(p,t.q,t.p))  return -1;
                else return -1;
            }
            else
            {
                if(Point.convex(t.p,q,p)) return 1;
                else return 0;
            }
        } 
}

class Point implements Comparable<Point>
{
    double x;
    double y;
    String Type;
    int index;
    
    public int compareTo(Point P) 
    {
        if(this.y > P.y) return -1;
        else if (this.y < P.y) return +1;
        else if (this.x < P.x) return -1;
        else if (this.x > P.x) return +1;
        else return  0; 
    }    
    
    static boolean below(Point p,Point q)
    {
        if(p.y < q.y) return true;
        else if(p.y == q.y)
        {
            if(p.x > q.x) return true;
        }
        return false;
    }
    
     static boolean IsConvex(Point p1, Point p2, Point p3)
    {
        double tmp;
        tmp = (p3.y-p1.y)*(p2.x-p1.x)-(p3.x-p1.x)*(p2.y-p1.y);
        if(tmp >0 ) return true;
        else return false;
    }
    
    static boolean convex(Point p,Point q,Point r)
    {
        double t;
        t= (p.x-q.x)*(p.y-r.y) - (p.x-r.x)*(p.y-q.y);
        if(t>0) return true;
        else return false;
    }
    
    static void Monotone(int N,Point point[],Vector<Vector<Point>> output,Vector<Vector<Point>> result)
    {
        Vector<Point> D=new Vector<>();
        //create priority queue of all the vertices 
        PriorityQueue<Point> PQ = new PriorityQueue<Point>();  
        for (int i = 0; i < N; i++) 
        {
            PQ.add(point[i]);
        }
        
        //find the type of all the vertices
        for(int i=1;i<N-1;i++)
        {
            if(below(point[i-1],point[i]) && below(point[i+1],point[i]))
            {
                if(convex(point[i-1],point[i],point[i+1]))
                    point[i].Type="start";
                else  point[i].Type="split";
            }
            else if(below(point[i],point[i-1]) && below(point[i],point[i+1]))
            {
                if(convex(point[i-1],point[i],point[i+1]))
                    point[i].Type="end";
                else  point[i].Type="marge";
            }
            else  point[i].Type="regular";
        }
        //find the type of the first vertex
        point[0].Type="start";
        //find the type of the last vertex
        if(below(point[N-2],point[N-1]) && below(point[0],point[N-1]))
        {
            if(convex(point[N-2],point[N-1],point[0]))
                point[N-1].Type="start";
            else  point[N-1].Type="split";
        }
        else if(below(point[N-1],point[N-2])&& below(point[N-1],point[0])) 
        {
            if(convex(point[N-2],point[N-1],point[0]))
                point[N-1].Type="end";
            else  point[N-1].Type="marge";    
        }
        else  point[N-1].Type="regular";
        
        TreeSet<Edge> T = new TreeSet<Edge>();
        int helper[]=new int[N];
        //calling appropriate procedure for all the vertices
        for(int i=0;i<N;i++)
        {
            Point v=PQ.remove();
            int j;
            for(j=0;j<N;j++)
            {
                if(point[j]==v) break;
            }
            switch(v.Type)
            {
                case "start":
                    start(v,point,T,helper,j,N);
                    break;
                    
                case "end":
                    end(v,point,T,helper,j,N,D);
                    break;
                    
                case "split":
                    split(v,point,T,helper,j,N,D);
                    break;
                    
                case "marge":
                    marge(v,point,T,helper,j,N,D);
                    break;
                
                case "regular":
                    regular(v,point,T,helper,j,N,D);
                    break;
            }
            
        }
        for(int i=0;i<D.size();i++)
        {
            if(i%2==0) 
            {
                if(D.elementAt(i+1).x < D.elementAt(i).x)
                {
                    Point temp=D.elementAt(i);
                    D.remove(i);
                    D.add(i+1, temp);
                }
            }
        }
        System.out.println("Added Diagonal : ");
        System.out.println();
        for(int i=0;i<D.size();i=i+2) 
        {
            System.out.println(D.elementAt(i).x+" "+D.elementAt(i).y+"   "+D.elementAt(i+1).x+" "+D.elementAt(i+1).y);
        }
        find_monotone_polygon(point,D,N,output);
        System.out.println();
        System.out.println("Monotone Polygon : ");
        System.out.println();
        for(int i=0;i<output.size();i++) 
        {
            for(int j=0;j<output.get(i).size();j++)
            {
                System.out.println(output.get(i).elementAt(j).x+" "+output.get(i).elementAt(j).y);
            }
            System.out.println();
        }
        for(int i=0;i<output.size();i++)
            triangulate(output.elementAt(i),result);
   
        
        System.out.println("Triangulation : ");
        System.out.println(result.size() + " Triangles");
        System.out.println();
        for(int i=0;i<result.size();i++) 
        {
            for(int j=0;j<result.get(i).size();j++)
            {
                System.out.println(result.get(i).elementAt(j).x+" "+result.get(i).elementAt(j).y);
            }
            System.out.println();
        }
        System.out.println("----------------------------------------");
    }
    
    static void find_monotone_polygon(Point[] point,Vector<Point> D,int N,Vector<Vector<Point>> output)
    {
        Vector<Point> monotone=new Vector<>();
        if(D.size()==0)
        {
            for(int i=0;i<N;i++) monotone.add(point[i]);
            output.add(monotone);
            return;
        }
        Point start=new Point();
        Point end=new Point();
        int endindex=-1;
        Vector<Integer> save=new Vector<>();
        for(int i=0;i<N;i++)
        {
            monotone.add(point[i]);
            if(D.contains(point[i])) 
            {
                int index=D.indexOf(point[i]);
                start=D.elementAt(index);
                D.remove(index);
                save.add(start.index);
                if(index<D.size())
                {
                    end=D.elementAt(index);
                    D.remove(index);
                }
                else 
                {
                    end=D.elementAt(index-1);
                    D.remove(index-1);
                }
                monotone.add(end);
                i=end.index;
                save.add(i);
            }
            if(i==N-1 || i==endindex) 
            {
                output.add(monotone);
                monotone=new Vector<>();
                endindex=end.index;
                if(D.size()==0) break;
                i=start.index-1;
                if(save.contains(i+1)) 
                {
                    int d=save.indexOf(i+1);
                    save.remove(d);
                    save.remove(d);
                }
            }
        }
        
        while(!save.isEmpty())
        {
            int s=save.firstElement();
            save.remove(0);
            int e=save.firstElement();
            save.remove(0);
            for(int i=s;i<=e;i++)
            {
                monotone.add(point[i]);
            }
            output.add(monotone);
            monotone=new Vector<>();
        }
    }
    
    
    /*static void find_monotone_polygon(Point[] point,Vector<Point> D,int N,Vector<Vector<Point>> output)
    {
        Vector<Point> monotone=new Vector<>();
        if(D.size()==0)
        {
            for(int i=0;i<N;i++) monotone.add(point[i]);
            output.add(monotone);
            return;
        }
        int endindex=-1;
        int k=0;
        int check=0,flag=0,savestart=0,saveend=0,save=0;
        Point start=D.firstElement();
        Point end=new Point();
        D.remove(start);
        for(int i=0;i<N;i++)
        {
            System.out.println("point "+point[i].x+" "+point[i].y);
            monotone.add(point[i]);
            if(point[i]==start) 
            {
                check++;
                if(check==2) 
                {
                    flag=1;
                    savestart=k;
                    saveend=save;
                }
                if(D.size()==0) break;
                end=D.firstElement();
                monotone.add(end);
                D.remove(end);
                save=end.index;
                i=end.index;
                k=start.index-1;
                if(D.size()!=0) 
                {
                    start=D.firstElement();
                    D.remove(start);
                }
            }
            if(i==N-1 || i==endindex) 
            {
                check=0;
                output.add(monotone);
                monotone=new Vector<>();
                endindex=end.index;
                System.out.println("start "+start.x+" "+k);
                if(D.size()==0) break;
                i=k;
            }
        }
        System.out.println(flag+"        "+savestart+" "+saveend);
        System.out.println(start.index+"        "+end.index);
        if(flag==1)
        {
            for(int i=savestart+1;i<=saveend;i++)
            {
                monotone.add(point[i]);
            }
            output.add(monotone);
            monotone=new Vector<>();
        }
        for(int i=start.index;i<=end.index;i++)
        {
            monotone.add(point[i]);
        }
        output.add(monotone);
        
    }
    */ 
    static void start(Point v,Point[] point ,TreeSet<Edge> T,int[] helper,int i,int N)
    {
        //Insert ei in T and set helper(ei) to vi.
        Edge edge=new Edge();
        edge.p=v;
        if(v.index==N-1) edge.q=point[0];
        else edge.q=point[v.index+1];
        T.add(edge);
        helper[i]=i;
    }
    
    static void end(Point v,Point point[],TreeSet<Edge> T,int[]helper,int i,int N,Vector<Point> D)
    {
        //if helper(ei−1) is a merge vertex
        if(point[helper[i-1]].Type=="marge") 
        {
            //then Insert the diagonal connecting vi to helper(ei−1) in D.
            D.add(v);
            D.add(point[helper[i-1]]);
            System.out.println("diagonal : ("+v.x+" "+v.y+")  ("+point[helper[i-1]].x+" "+point[helper[i-1]].y+")");
        }
        //Delete ei−1 from T.
        Edge edge=new Edge();
        if(v.index==0) edge.q=point[N-1]; 
        else edge.p=point[v.index-1];
        edge.q=v;
        T.remove(edge);
    }
    
    static void split(Point v,Point[] point,TreeSet<Edge> T,int[]helper,int i,int N,Vector<Point> D)
    {
       // Search in T to find the edge ej directly left of vi.
        /*Iterator<Edge> iterator = T.iterator();
        Edge e1=new Edge();
        double max=-999;
        int ej=0;
        while(iterator.hasNext())
        {
            e1=iterator.next();
            if(e1.p.x < v.x && e1.p.x > max )
            {
                max=e1.p.x;
                ej=e1.p.index;
                left=helper[e1.p.index];
            }
        }*/
        Edge edg=new Edge();
        edg.p=v;
        edg.q=v;
        edg=T.lower(edg);
        int left=helper[edg.p.index];
        
        // Insert the diagonal connecting vi to helper(ej) in D.
        D.add(v);
        D.add(point[left]);
       // helper(ej)=vi
        helper[edg.p.index]=i;
       // Insert ei in T and set helper(ei) to vi.
        Edge edge=new Edge();
        edge.p=v;
        if(v.index==N-1) edge.q=point[0]; 
        else edge.q=point[v.index+1];
        T.add(edge);
        helper[i]=i;
    }
    
    static void marge(Point v,Point[] point,TreeSet<Edge> T,int[]helper,int i,int N,Vector<Point> D)
    {
        //if helper(ei−1) is a merge vertex
        if(point[helper[i-1]].Type=="marge") 
        {
            //then Insert the diagonal connecting vi to helper(ei−1) in D.
            D.add(v);
            D.add(point[helper[i-1]]);
        }
        //Delete ei−1 from T.
        Edge edge=new Edge();
        if(v.index==0) edge.p=point[N-1]; 
        else edge.p=point[v.index-1];
        edge.q=v;
        T.remove(edge);
      
        Edge edg=new Edge();
        edg.p=v;
        edg.q=v;
        edg=T.lower(edg);
        int left=helper[edg.p.index];
        //if helper(ej) is a merge vertex then Insert the diagonal connecting vi to helper(ej) in D.
        if(point[left].Type=="marge") 
        {
            D.add(v);
            D.add(point[left]);
        }
        // helper(ej)=vi
        helper[edg.p.index]=i;
    }
    
    static void regular(Point v,Point[] point,TreeSet<Edge> T,int[]helper,int i,int N,Vector<Point> D)
    {
        // if the interior of P lies to the right of vi
        if(below(v,point[i-1]))
        {
            // then if helper(ei−1) is a merge vertex
            if(point[helper[i-1]].Type=="marge") 
            {
                // then Insert the diagonal connecting vi to helper(ei−1) in D.
                D.add(v);
                D.add(point[helper[i-1]]);
            }
            //Delete ei−1 from T.
            Edge edge=new Edge();
            if(v.index==0) edge.q=point[N-1]; 
            else edge.p=point[v.index-1];
            edge.q=v;
            T.remove(edge);
            // Insert ei in T and set helper(ei) to vi.
            edge.p=v;
            if(v.index==N-1) edge.q=point[0]; 
            else edge.q=point[v.index+1];
            T.add(edge);
            helper[i]=i;
        }
        // else Search in T to find the edge ej directly left of vi.
        else 
        {
            Edge edg=new Edge();
            edg.p=v;
            edg.q=v;
            edg=T.lower(edg);
            int left=helper[edg.p.index];
            
            //if helper(ej) is a merge vertex
            if(point[left].Type=="marge") 
            {
                // then Insert the diagonal connecting vi to helper(ej) in D.
                D.add(v);
                D.add(point[left]);
            }
            // helper(ej)=vi
            helper[edg.p.index]=i;
        }
    }
    
    static void triangulate(Vector<Point> Monotone, Vector<Vector<Point>> triangles)
    {
        Vector<Point> triangle=new Vector<>();

        int numpoints = Monotone.size();
        
        if(numpoints == 3)
        {
            triangles.add(Monotone);
            return ;
        }

        int topindex = 0;  
        int bottomindex=0;
        int i=0;
        for(i=1;i<numpoints;i++)
        {
            if(below(Monotone.elementAt(i), Monotone.elementAt(bottomindex))) bottomindex = i;
            if(below(Monotone.elementAt(topindex),Monotone.elementAt(i))) topindex = i;
        }
        int vertextypes[] = new int[numpoints];
        int priority[] = new int[numpoints];

        //Merge the vertices on the left chain and the vertices on the right chain of P
        priority[0] = topindex;
        vertextypes[topindex] = 0;
        int leftindex = topindex+1;
        if(leftindex>=numpoints) leftindex = 0;
        int rightindex = topindex-1;
        if(rightindex<0) rightindex = numpoints-1;
        for(i=1;i<(numpoints-1);i++)
        {
           if(leftindex==bottomindex)
            {
                priority[i] = rightindex;
                rightindex--;
                if(rightindex<0) rightindex = numpoints-1;
                vertextypes[priority[i]] = -1;
            }
            else if(rightindex==bottomindex)
            {
                priority[i] = leftindex;
                leftindex++;
                if(leftindex>=numpoints) leftindex = 0;
                vertextypes[priority[i]] = 1;
            }
            else
            {
                if(below(Monotone.elementAt(leftindex),Monotone.elementAt(rightindex)))
                {
                    priority[i] = rightindex;
                    rightindex--;
                    if(rightindex<0) rightindex = numpoints-1;
                    vertextypes[priority[i]] = -1;
                }
                else
                {
                    priority[i] = leftindex;
                    leftindex++;
                    if(leftindex>=numpoints) leftindex = 0;
                    vertextypes[priority[i]] = 1;
                }
            }
        }
        priority[i] = bottomindex;
        vertextypes[bottomindex] = 0;

        //Initialize an empty stack S 
        int S[] = new int[numpoints];
        int stackptr = 0;
        //push u1 and u2 onto S.
        S[0] = priority[0];
        S[1] = priority[1];

        int vindex=0;
        stackptr = 2;
        //for j=3 to n−1
        for(i=2;i<(numpoints-1);i++)
        {
            vindex = priority[i];
            //if uj and the vertex on top of S are on different chains
            if(vertextypes[vindex]!=vertextypes[S[stackptr-1]])
            {
                for(int j=0;j<(stackptr-1);j++)
                {
                    if(vertextypes[vindex]==1)
                    {
                        triangle.add(Monotone.elementAt(S[j+1]));
                        triangle.add(Monotone.elementAt(S[j]));
                        triangle.add(Monotone.elementAt(vindex));
                    }
                    else
                    {
                        triangle.add(Monotone.elementAt(S[j]));
                        triangle.add(Monotone.elementAt(S[j+1]));
                        triangle.add(Monotone.elementAt(vindex));
                    }
                    triangles.add(triangle);
                    triangle=new Vector<>();
                }
                //Push uj−1 and uj onto S.
                S[0] = priority[i-1];
                S[1] = priority[i];
                stackptr = 2;
            }
            else
            {
                stackptr--;
                while(stackptr>0)
                {
                    if(vertextypes[vindex]==1)
                    {
                        if(convex(Monotone.elementAt(S[stackptr-1]),Monotone.elementAt(S[stackptr]),Monotone.elementAt(vindex)))
                        {
                            triangle.add(Monotone.elementAt(vindex));
                            triangle.add(Monotone.elementAt(S[stackptr-1]));
                            triangle.add(Monotone.elementAt(S[stackptr]));
                            triangles.add(triangle);
                            triangle=new Vector<>();
                            stackptr--;
                        }
                        else
                            break;
                    }
                    else
                    {
                        if(convex(Monotone.elementAt(S[stackptr]),Monotone.elementAt(S[stackptr-1]),Monotone.elementAt(vindex)))
                        {
                            triangle.add(Monotone.elementAt(vindex));
                            triangle.add(Monotone.elementAt(S[stackptr]));
                            triangle.add(Monotone.elementAt(S[stackptr-1]));
                            triangles.add(triangle);
                            triangle=new Vector<>();
                            stackptr--;
                        }
                        else
                            break;
                    }
                }
                stackptr++;
                //Push uj onto S.
                S[stackptr] = vindex;
                stackptr++;
            }
        }
        //Add diagonals from un to all stack vertices except the first and the last one.
        vindex = priority[i];
        for(int j=0;j<(stackptr-1);j++)
        {
            if(vertextypes[S[j+1]]==1)
            {
                triangle.add(Monotone.elementAt(S[j]));
                triangle.add(Monotone.elementAt(S[j+1]));
                triangle.add(Monotone.elementAt(vindex));
            }
            else
            {
                triangle.add(Monotone.elementAt(S[j+1]));
                triangle.add(Monotone.elementAt(S[j]));
                triangle.add(Monotone.elementAt(vindex));
            }
            triangles.add(triangle);
            triangle=new Vector<>();
        }
        
        return ;
    }
}   
    
public class Offline_3 {

    /**
     * @param args the command line arguments
     */
    public static Scanner in;
    public static void main(String[] args) throws FileNotFoundException {
        in = new Scanner(new File("input.txt"));
        int T=in.nextInt();
        for(int l=0;l<T;l++)
        {
            int N=in.nextInt();
            Point point[]=new Point[N];
            for (int i=0;i<N;i++)
            {
                point[i]=new Point();
                point[i].x=in.nextInt();
                point[i].y=in.nextInt();
                point[i].index=i;
            }
            Vector<Vector<Point>> output = new Vector<Vector<Point>>();
            Vector<Vector<Point>> result = new Vector<Vector<Point>>();
            Point.Monotone(N,point,output,result);
            System.out.println();
            System.out.println();
        }
    }  
}
