//26/05/2020
//Tortoise Name: Golam Rabbi Nazum
/*
->Sometimes being a stubborn is helpful!
->Confidence is the main weapon!
->Happy coding!
*/
#include<bits/stdc++.h>
#include<cmath>
#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;

#define debug printf("Debug\n");
#define en '\n'
#define all(x) x.begin() , x.end()

#define ll long long
#define ull unsigned long long
#define dob double

#define pi acos(-1.0)
#define mod 1000000007

#define _fastio  ios_base:: sync_with_stdio(false); cin.tie(0); cout.tie(0);
#define in_fre  freopen("in.txt", "r", stdin);
#define out_fre freopen("out.txt", "w", stdout);

struct Node{
    ll x,y;
};
struct cmp{
    bool operator()(const Node &a,const Node &b)const{
        if(a.x==b.x) return a.y<b.y;
        else return a.x>b.x;
    }
};
struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        // http://xorshift.di.unimi.it/splitmix64.c
        //only submit on GNU C++20 (64)
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};
//STL Algos
template<typename T> using Pbds=tree<T,null_type,less<ll>,rb_tree_tag,tree_order_statistics_node_update>;

//Number Theory
class NumberTheory{
public:
    ll n;
    ll mindiv[200000],fact[1000006];
    bool chk[99999995];
    ll dp[5005][5005];
    
    ll GetMod(ll x,ll m);
    ll BigMod (ll b,ll p,ll m){if (p == 0) return 1;if (p%2 == 0){ll s = BigMod(b,p/2,m);return ((s%m)*(s%m))%m;}return ((b%m)*(BigMod(b,p-1,m)%m))%m;}
    ll ModInv (ll b,ll m){return BigMod(b,m-2,m);}
    ll gcd(ll a,ll b){if(a<0)return gcd(-a,b);if(b<0)return gcd(a,-b);return (b==0)?a:gcd(b,a%b);}
    ll lcm(ll a,ll b) {if(a<0)return lcm(-a,b);if(b<0)return lcm(a,-b);return a*(b/gcd(a,b));}
    void sieve(ll x);
    bool IsPrime(ll x);
    bool IsPrime1(ll x){return !chk[x];}
    void FactSieve(ll x);
    vector<ll> getFactorization(ll x);
    ll comb(ll n,ll r);
    void gen_fact(ll n){fact[0]=1;for(ll i=1;i<=n;i++){fact[i]=fact[i-1]*i;fact[i]%=mod;}}
    ll get_comb(ll n,ll r){ll res=fact[n];res*=ModInv(fact[r],mod);res%=mod;res*=ModInv(fact[n-r],mod);res%=mod;return res;}
};
//Geometry
class Geometry{
public:
struct CartesianPoint{dob x,y;};
struct PolarPoint{dob r,ang;};
    dob SQR(dob x);
    dob DOT(dob x1,dob y1,dob x2,dob y2);
    dob CROSS(dob x1,dob y1,dob x2,dob y2);
    dob PolarR(dob x,dob y);
    dob PolarAngle(dob x,dob y);
    dob CartesianX(dob r,dob ang);
    dob CartesianY(dob r,dob ang);
    dob CartesianDis(dob x1,dob y1,dob x2,dob y2);
    dob PolarDis(dob r1,dob ang1,dob r2,dob ang2);
    dob MiddlePoint(dob x1,dob x2);
    dob InternalSection(dob x1,dob x2,dob m1,dob m2);
    dob ExternalSection(dob x1,dob x2,dob m1,dob m2);
    dob AreaOfTriangle(dob x1,dob y1,dob x2,dob y2,dob x3,dob y3);
    dob GradAngle(dob x1,dob y1,dob x2,dob y2);
    dob Grad(dob x1,dob y1,dob x2,dob y2);
    dob Line(dob m,dob c,dob x,dob y);
    dob Line(dob m,dob x1,dob y1,dob x,dob y);
    dob Line(dob x1,dob y1,dob x2,dob y2,dob x,dob y);
    dob OrthogonalLine(dob x1,dob y1,dob x2,dob y2,dob lx,dob ly,dob ox,dob oy);
    dob LineOrthogonalLineIntersecX(dob x1,dob y1,dob x2,dob y2,dob ox,dob oy);
    dob LineOrthogonalLineIntersectY(dob x1,dob y1,dob x2,dob y2,dob ox,dob oy);
    bool IsParallel(dob x1,dob y1,dob x2,dob y2);
    bool IsOrthogonal(dob x1,dob y1,dob x2,dob y2);
    dob IntersecX(dob a1,dob b1,dob c1,dob a2,dob b2,dob c2);
    dob IntersecY(dob a1,dob b1,dob c1,dob a2,dob b2,dob c2);
    dob PointToLineDis(dob a,dob b,dob c,dob x,dob y);
    dob PointToLineDis(dob x1,dob y1,dob x2,dob y2,dob x,dob y);
    dob LineToLineDis(dob a,dob b,dob c1,dob c2);
    ll isLeft(CartesianPoint a,CartesianPoint b,CartesianPoint c);

};

//Bitwise operations
class Bitwise_operation{
public:
    ull a[65];
    ull A,B;
    bool check(ull N,ull pos){ return (bool)(N & (1LL<<pos)); }
    ull Set1(ull N,ull pos){ return N=N | (1LL<<pos); }
    ull Set0(ull N,ull pos){ return N=N & ~(1LL<<pos); }
    ull flip(ull N,ull pos){ return N=N ^ (1LL<<pos); }
    bool isPowerOf2(ull n){return !(n&(n-1))?true:false;}
    void swap(){ A = A ^ B;B = A ^ B;A = A ^ B;}
    ull FindLargestPowerof2(ull N){N=N|N>>1;N=N|N>>2;N=N|N>>4;N=N|N>>8;N=N|N>>16;N=N|N>>32; return(N+1)>>1;}
    ull GetRightmost1(ull x){return (x & (-x));}
    ull GetRightmost1spos(ull x){return __builtin_ffs(x)-1;}
    ll CountNumberof1s(ull x){ll cnt=0;while( x ){x = x&(x-1);cnt++;}return cnt;}
    ll CountLeadingZeros(ull x){return __builtin_clz(x);}
    ll CountTrailingZeros(ull x){return __builtin_ctz(x);}
    void possibleSubsets(ll N);

};

//Binary Search
class Binary_Search{
public:
    ll ar[200005];
    ll Manual_upper(ll l,ll r,ll val);
    ll Manual_lower(ll l,ll r,ll val);
};

//Graph Theory
class GraphTheory{
public:
    struct Edge {
        ll a, b, cost;
    };
    vector<Edge>allEdges;
    ll nodes,edges,INF=1e18;
    vector<ll>graph[200005],weight[200005];
    vector<ll>d,p;
    vector<bool> u;
    ll vis[200005];ll cnt[200005];ll par[200005],lvl[200005],parse_table[200005][20];
    void CLR(ll n);
    bool isLeaf(ll node);
    void Connect(ll u,ll v);
    void weighted_connect(ll u,ll v,ll w);
    void bfs(ll src);
    void dfs(ll src);
    void CountSubtreeSize(ll cur);
    ll dirx[5]={1,-1,0,0};
    ll diry[5]={0,0,1,-1};
    void GirdBfs(ll i,ll j);
    void GridDfs(ll i,ll j);
    void generate_parse_table();
    ll lca(ll p,ll q);
    void dijkstra(ll s);
    void dijkstra1(ll s);
    void BellManFord(ll s);
    vector<ll> FindNegativeCycle(ll s);
    bool SPFA(ll s); 
    vector<ll> restore_path(ll s, ll t);
    
    vector<bool> visited;
    vector<ll> tin, low;
    vector<pair<ll,ll>>bridges;
    ll timer = 0;
    void find_bridges(ll n);
    void dfs_time(ll v, ll p);
    
    vector<ll>articulation_points;
    void find_Articulation_Point(ll n);
    void dfs_time1(ll v, ll p);
};
//Data Structure
class DataStructure_DSU{
public:
    ll parent[200005];ll Size[200005];ll Rank[200005];
    void CLR(ll n);
    void make_set(ll v);
    ll find_set(ll v);
    void union_sets(ll a, ll b);
    ll Path_Compressed_find_set(ll v);
    void union_sets_considering_size(ll a, ll b);
    void union_sets_considering_rank(ll a, ll b);
    void union_sets_considering_size_without_path_compression(ll a, ll b);
};

class DataStructure_Trie{
public:
    bool f;
    ll tot;
    ll tri[400005][30],emark[400005],cnt[400005],depth[400005];
    string s;
    void CLR(ll n);
    void build();
    ll Find();
};
class DataStructure_Trie1{
public:
    struct node{
        ll cnt;
        node* nxt[26];
    
        node(){
            cnt=0;
            for(int i=0;i<26;i++){
                nxt[i]=nullptr;
            }
        }
    };
    string s;
    node* root=new node();
    void build();
    ll Find();
};

class DataStructure_Segment_Tree{
public:
    ll n,tot;
    ll el[100005],segtree[400005];
    void CLR();
    void build(ll node,ll b,ll e);
    void up(ll node,ll b,ll e,ll p);
    void up1(ll node,ll b,ll e,ll l,ll r);
    void query(ll node,ll b,ll e,ll p);
    void query1(ll node,ll b,ll e,ll l,ll r);
};

class DataStructure_SQRT_DCOM{
public:

};

class Hashing{
public:
    string s;
    ll n;
    char a[500005];
    ull forwdhash[500005],revhash[500005];
    ull po[500005],revpo[500005];
    void build();
    bool isPlain(ll l,ll r);
};

class Miscellaneous{
public:
    ll lis(vector<ll> const& a);//nlog(n)
    bool issubsequence(string& s1, string& s2);//O(max(l1,l2))
};

/*Containers*/
priority_queue<ll,deque<ll>,greater<ll>>pq;
set<Node,cmp>st1;
map<Node,ll,cmp>mp1;
priority_queue<Node,deque<Node>,cmp>pq1;
unordered_map<ll,ll,custom_hash>ump1;
unordered_set<ll,ll,custom_hash>ust1;
Pbds<ll>ost;





int main(){
    _fastio
    ll cas;
    cin>>cas;
    while(cas--){
        
    }


    return 0;
}

//NumberTheory
void NumberTheory::sieve(ll x){
    for(ll i=2;i*i<=x;){
        if(chk[i]==0){
            for(ll j=i*i;j<=x;j+=i){
                chk[j]=1;
            }
        }
        if(i==2) i++;
        else i+=2;
    }
}
bool NumberTheory::IsPrime(ll x){
    bool f=true;
    for(ll i=2;i*i<=x;i++){
        if(!x%i){
            f=false;break;
        }
    }
    return f;
}
ll GetMod(ll x,ll m){
    if(x>=0) return x%m;
    else{
        ll y=abs(x)/m;x+=(y*m+m+m);
        return x%m;
    }
}
void NumberTheory::FactSieve(ll x){
    mindiv[1] = 1;
    for (ll i=2; i<x; i++) mindiv[i] = i;
    for (ll i=4; i<x; i+=2) mindiv[i] = 2;
    for (ll i=3; i*i<x; i++){
        if (mindiv[i] == i){
            for (ll j=i*i; j<x; j+=i)
                if (mindiv[j]==j)
                    mindiv[j] = i;
        }
    }
}
vector<ll> NumberTheory::getFactorization(ll x){
    vector<ll> ret;
    while (x != 1){
        ret.emplace_back(mindiv[x]);
        x = x / mindiv[x];
    }
    return ret;
}

ll NumberTheory::comb(ll n,ll r){
    if(n<0 or r<0) return 0;
    if(dp[n][r]!=-1) return dp[n][r];
    if(n<r) return dp[n][r]=0;
    if(r == 0) return dp[n][r]=1;
    if(r == 1) return dp[n][r]=n;
    if(n == 1)return dp[n][r]=1;
    return dp[n][r]=comb(n-1,r-1)+comb(n-1,r);
}

//Bitwise operations
void Bitwise_operation::possibleSubsets(ll N){
    for(ll i = 0;i <(1 << N); ++i){
        for(ll j = 0;j < N;++j)
            if(check(i,j))
                cout << a[j]<<' ';
        cout << endl;
    }
}

//BinarySearch
ll Binary_Search::Manual_upper(ll l,ll r,ll val){
    ll x=val;
    ll pos=0;
    while(1){
        ll mid=(l+r)>>1;
        if(r-l<=1){
            if(ar[r]==x) pos=r;
            else if(ar[l]==x) pos=l;
            else pos=l;
            break;
        }
        else if(x>=ar[mid]) l=mid;
        else r=mid;
    }
    return pos;
}
ll Binary_Search::Manual_lower(ll l,ll r,ll val){
    ll x=val;
    ll pos=0;
    while(1){
        ll mid=(l+r)>>1;
        if(r-l<=1){
            if(ar[l]==x) pos=l;
            else if(ar[r]==x) pos=r;
            else pos=l;
            break;
        }
        else if(x<=ar[mid]) r=mid;
        else l=mid;
    }
    return pos;
}

//Graph Theory
void GraphTheory::CLR(ll n){
    for(ll i=0;i<=n+1;i++){
        graph[i].clear();
        weight[i].clear();
        vis[i]=0;
    }
    nodes=n;
    d.assign(n+1, 1e18);
    p.assign(n+1, -1);
    u.assign(n+1, false);
    INF=1e18;
}
bool GraphTheory::isLeaf(ll node){
    if(graph[node].size()==1) return true;
    else return false;
}

void GraphTheory::Connect(ll u,ll v){
    graph[u].push_back(v);
    graph[v].push_back(u);
}
void GraphTheory::bfs(ll src){
    vis[src]=1;
    deque<ll>dq;
    dq.push_back(src);
    while(!dq.empty()){
        src=dq.front();
        dq.pop_front();
        for(ll i=0;i<graph[src].size();i++){
            ll adj=graph[src][i];
            if(vis[adj]==0){
                vis[adj]=1;
                dq.push_back(adj);
            }
        }
    }
}
void GraphTheory::dfs(ll src){
    vis[src]=1;
    for(ll i=0;i<graph[src].size();i++){
        ll adj=graph[src][i];
        if(!vis[adj]){
            dfs(adj);
        }
    }
}
void GraphTheory::CountSubtreeSize(ll cur){
    vis[cur]=1;
    cnt[cur]++;
    for(ll i=0;i<graph[cur].size();i++){
        ll adj=graph[cur][i];
        if(!vis[adj]){
            CountSubtreeSize(adj);
            cnt[cur]+=cnt[adj];
        }
    }
}

void GraphTheory::generate_parse_table(){
    for(ll i=1;i<=nodes;i++){
        parse_table[i][0]=par[i];
    }
    for(ll j=1;j<=20;j++){
        for(ll i=1;i<=nodes;i++){
            if(parse_table[i][j-1]!=-1)
                parse_table[i][j]=parse_table[parse_table[i][j-1]][j-1];
            else
                parse_table[i][j]=-1;
        }
    }
}

ll GraphTheory::lca(ll p,ll q){
   for(ll i=20;i>=0;i--){
        if(lvl[p]==lvl[q]) break;
        if(lvl[p]>lvl[q]){
            if(lvl[p]-(1<<i)>=lvl[q]) p=parse_table[p][i];
        }
        else{
            if(lvl[q]-(1<<i)>=lvl[p]) q=parse_table[q][i];
        }
    }

    if(p==q) return p;

    for(ll i=20;i>=0;i--){
        if(parse_table[p][i]!=-1&&parse_table[p][i]!=parse_table[q][i]){
            p=parse_table[p][i];
            q=parse_table[q][i];
        }
    }
    return parse_table[p][0];
}

void GraphTheory::weighted_connect(ll u,ll v,ll w){
    graph[u].push_back(v);
    graph[v].push_back(u);
    weight[u].push_back(w);
    weight[v].push_back(w);
}

void GraphTheory::dijkstra(ll s){
    //n^2+m
    d[s] = 0;
    for (ll i = 0; i < nodes; i++) {
        ll v = -1;
        for (ll j = 0; j < nodes; j++) {
            if (!u[j] && (v == -1 || d[j] < d[v]))
                v = j;
        }

        if (d[v] == INF)
            break;
        u[v] = true;
        for (ll j=0;j<graph[v].size();j++) {
            ll to = graph[v][j];
            ll w = weight[v][j];

            if (d[v] + w < d[to]) {
                d[to] = d[v] + w;
                p[to] = v;
            }
        }
    }
}

void GraphTheory::dijkstra1(ll s) {
    //nlong+m
    d[s] = 0;
    using pii = pair<ll, ll>;
    priority_queue<pii, vector<pii>, greater<pii>> q;
    q.push({0, s});
    while (!q.empty()) {
        ll v = q.top().second;
        ll d_v = q.top().first;
        q.pop();
        if (d_v != d[v])
            continue;
        for (ll j=0;j<graph[v].size();j++) {
            ll to = graph[v][j];
            ll w = weight[v][j];
            if (d[v] + w < d[to]) {
                d[to] = d[v] + w;
                p[to] = v;
                q.push({d[to], to});
            }
        }
    }
}

void GraphTheory::BellManFord(ll s){
    //n*m
    //works for negative weight edge
    //no negative cycle
    d[s] = 0;
    for (;;) {
        bool any = false;
        for (auto e : allEdges)
            if (d[e.a] < INF)
                if (d[e.b] > d[e.a] + e.cost) {
                    d[e.b] = d[e.a] + e.cost;
                    p[e.b] = e.a;
                    any = true;
                }
        if (!any)
            break;
    }
}

vector<ll> GraphTheory::FindNegativeCycle(ll s){
    //n*m
    //works for negative weight edge
    //negative cycle detection
    d[s] = 0;
    ll x;
    for (ll i = 0; i < nodes; ++i) {
        x = -1;
        for (auto e : allEdges)
            if (d[e.a] < INF)
                if (d[e.b] > d[e.a] + e.cost) {
                    d[e.b] = max(-INF, d[e.a] + e.cost);
                    p[e.b] = e.a;
                    x = e.b;
                }
    }

    vector<ll> path;
    if (x == -1){
        cout << "No negative cycle from " << s;
        return path;
    }
    else {
        ll y = x;
        for (ll i = 0; i < nodes; ++i)
            y = p[y];
        for (ll cur = y;; cur = p[cur]) {
            path.push_back(cur);
            if (cur == y && path.size() > 1)
                break;
        }
        reverse(path.begin(), path.end());

        cout << "Negative cycle: ";
        return path;
    }
}

bool GraphTheory::SPFA(ll s){
    //shortest path faster algorithm
    //n*m
    //if there is negative cycle that can be detected otherwise single source shortest path
    vector<ll> cnt(nodes+1, 0);
    vector<bool> inqueue(nodes+1, false);
    queue<ll> q;

    d[s] = 0;
    q.push(s);
    inqueue[s] = true;
    while (!q.empty()) {
        ll v = q.front();
        q.pop();
        inqueue[v] = false;

        for (ll i=0;i<graph[v].size();i++) {
            ll to = graph[v][i];
            ll w = weight[v][i];
            if (d[v] + w < d[to]) {
                d[to] = d[v] + w;
                if (!inqueue[to]) {
                    q.push(to);
                    inqueue[to] = true;
                    cnt[to]++;
                    if (cnt[to] > nodes)
                        return false;  // negative cycle
                }
            }
        }
    }
    return true;

}



vector<ll> GraphTheory::restore_path(ll s, ll t) {
    vector<ll> path;
    if(d[t]==INF) return path;

    for (ll v = t; v != s; v = p[v])
        path.push_back(v);
    path.push_back(s);

    reverse(path.begin(), path.end());
    return path;
}



void GraphTheory::find_bridges(ll n) {
    //O(n+m) finding all bridges together
    timer = 0;
    visited.assign(n, false);
    tin.assign(n, -1);
    low.assign(n, -1);
    for (ll i = 0; i < n; ++i) {
        if (!visited[i])
            dfs_time(i,-1);
    }
}

void GraphTheory::dfs_time(ll v, ll p) {
    visited[v] = true;
    tin[v] = low[v] = timer++;
    for (ll to : graph[v]) {
        if (to == p) continue;
        if (visited[to]) {
            low[v] = min(low[v], tin[to]);
        } else {
            dfs_time(to, v);
            low[v] = min(low[v], low[to]);
            if (low[to] > tin[v])
                bridges.emplace_back(make_pair(v,to));
        }
    }
}

void GraphTheory::find_Articulation_Point(ll n) {
    //O(n+m)
    timer = 0;
    visited.assign(n, false);
    tin.assign(n, -1);
    low.assign(n, -1);
    for (ll i = 0; i < n; ++i) {
        if (!visited[i])
            dfs_time1 (i,-1);
    }
}

void GraphTheory::dfs_time1(ll v, ll p) {
    visited[v] = true;
    tin[v] = low[v] = timer++;
    ll children=0;
    for (ll to : graph[v]) {
        if (to == p) continue;
        if (visited[to]) {
            low[v] = min(low[v], tin[to]);
        } else {
            dfs_time1(to, v);
            low[v] = min(low[v], low[to]);
            if (low[to] >= tin[v] && p!=-1)
                articulation_points.emplace_back(v);
            ++children;
        }
    }
    if(p == -1 && children > 1)
        articulation_points.emplace_back(v);
}





//Data Structure
/*
Using DSU.
1.we can find number of connected components of a graph
2.we can check either two vertex are in same component or different
3.we can find how many vertex is in each component
4.we can find minimum spanning tree
5.using find_parent function without path compression we can find any ancestor of a vertex. that means we can find lowest common ancestor
6.Also, it's worth mentioning that DSU with union by size / rank, but without path compression works in O(logn) time per query.
7.
*/
void DataStructure_DSU::CLR(ll n){
    for(ll i=0;i<=n+1;i++){
        parent[i]=i;
        Rank[i]=0;
        Size[i]=1;
    }
}
void DataStructure_DSU::make_set(ll v) {
    parent[v] = v;
}
ll DataStructure_DSU::find_set(ll v) {
    if (v == parent[v])
        return v;
    return find_set(parent[v]);
}
void DataStructure_DSU::union_sets(ll a, ll b) {
    a = find_set(a);
    b = find_set(b);
    if (a != b)
        parent[b] = a;
}
ll DataStructure_DSU::Path_Compressed_find_set(ll v) {
    if (v == parent[v])
        return v;
    return parent[v]=Path_Compressed_find_set(parent[v]);
}
void DataStructure_DSU::union_sets_considering_size(ll a, ll b) {
    a = Path_Compressed_find_set(a);
    b = Path_Compressed_find_set(b);
    if (a != b) {
        if (Size[a] < Size[b])
            swap(a, b);
        parent[b] = a;
        Size[a] += Size[b];
    }
}
void DataStructure_DSU::union_sets_considering_rank(ll a, ll b){
    a = Path_Compressed_find_set(a);
    b = Path_Compressed_find_set(b);
    if (a != b) {
        if (Rank[a] < Rank[b])
            swap(a, b);
        parent[b] = a;
        if (Rank[a] == Rank[b])
            Rank[a]++;
    }
}
void DataStructure_DSU::union_sets_considering_size_without_path_compression(ll a, ll b){
    //log(n)
    a = find_set(a);
    b = find_set(b);
    if (a != b) {
        if (Size[a] < Size[b])
            swap(a, b);
        parent[b] = a;
        Size[a] += Size[b];
    }
}




/*
*We can check whether an string is a prefix/suffix of any set of strings OR any string from a set of strings is a prefix/suffix of that string or not.
*We can find how many times a string appears as a prefix/suffix in a set of strings.
*We can find longest/shortest common prefix/suffix
*We can find longest common substring
*/

void DataStructure_Trie::CLR(ll n){
    for(ll i=0;i<=n+1;i++){
        memset(tri[i],0,sizeof(tri[i]));
        emark[i]=0;
        cnt[i]=0;
    }
    s.clear();
    tot=0;
}

void DataStructure_Trie::build(){
    ll now=0;
    ll len=s.length();
    for(ll i=0;i<len;i++){
        ll x=s[i]-'0';
        if(tri[now][x]==0){
            tot++;
            tri[now][x]=tot;
            depth[tot]=depth[now]+1;
            now=tot;
        }
        else{
            now=tri[now][x];
        }
        cnt[now]++;
        if(i==len-1) emark[now]++;
    }
}

ll DataStructure_Trie::Find(){
    ll ans=0;
    ll now=0;
    ll len=s.length();
    for(ll i=0;i<len;i++){
        ll x=s[i]-'0';
        if(tri[now][x]==0){
            return 0;
        }
        else{
            now=tri[now][x];
        }
        if(i==len-1){
            ans=cnt[now];
        }
    }
    return ans;
}
void DataStructure_Trie1::build(){
    node* now=root;
    //now->cnt++;
    for(ll i=0;i<s.length();i++){
        ll x=s[i]-'a';
        if(!now->nxt[x]){
            now->nxt[x]=new node();
        }
        now=now->nxt[x];
        now->cnt++;
    }
}

ll DataStructure_Trie1::Find(){
    ll res=0;
    node* now=root;
    for(auto c:s){
        ll x=c-'a';
        if(!now->nxt[x]) break;
        now=now->nxt[x];
        res+=2*(now->cnt);
    }
    
    return res;
}
/*
A lot of works we can do . :p
*/
void DataStructure_Segment_Tree::CLR(){
    memset(segtree,0,sizeof(segtree));
}
void DataStructure_Segment_Tree::build(ll node,ll b,ll e){
    if(b==e){
        segtree[node]=el[b];
        return;
    }
    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;
    build(left,b,mid);
    build(right,mid+1,e);
    segtree[node]=segtree[left]+segtree[right];
    return;

}

void DataStructure_Segment_Tree::up(ll node,ll b,ll e,ll p){
     if(b==e){
        if(b==p)
            segtree[node]=el[b];
        return;
    }
    if(b>p||e<p) return;

    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;
    up(left,b,mid,p);
    up(right,mid+1,e,p);
    segtree[node]=segtree[left]+segtree[right];
    return;
}
void DataStructure_Segment_Tree::up1(ll node,ll b,ll e,ll l,ll r){
    if(b==e){
        if(b>=l&&e<=r)
            segtree[node]++;
        return;
    }
    if(b>r||e<l) return;

    if(b>=l&&e<=r){
        segtree[node]++;
        return;
    }

    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;
    up1(left,b,mid,l,r);
    up1(right,mid+1,e,l,r);
}

void DataStructure_Segment_Tree::query(ll node,ll b,ll e,ll p){
    if(b==e){
        if(b==p)
            tot+=segtree[node];
        return;
    }
    if(b>p||e<p) return;

    tot+=segtree[node];

    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;
    query(left,b,mid,p);
    query(right,mid+1,e,p);
    return;
}
//void query(ll node,ll b,ll e,ll l,ll r);

/*
struct Node{
    ll x,y,val;
    bool f;
};
void build1(ll node,ll b,ll e){
    if(b==e){
        segtree[node].val=a[b];
        segtree[node].f=true;
        segtree[node].x=0;
        return;
    }
    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;
    segtree[node].val=0;
    segtree[node].x=0;
    segtree[node].f=false;
    build1(left,b,mid);
    build1(right,mid+1,e);
    
    return;

}
void up1(ll node,ll b,ll e,ll l,ll r,ll va){
    if(b==e){
        if(b>=l&&e<=r)
            segtree[node].x+=va;
            segtree[node].f=true;
        return;
    }
    if(b>r||e<l) return;

    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;

    if(segtree[node].f){
        segtree[left].f=true;
        segtree[left].x+=segtree[node].x;
        segtree[right].f=true;
        segtree[right].x+=segtree[node].x;

        segtree[node].x=0;
        segtree[node].f=false;
    }

    if(b>=l&&e<=r){
        segtree[node].x+=va;
        segtree[node].f=true;
        return;
    }

    
    up1(left,b,mid,l,r,va);
    up1(right,mid+1,e,l,r,va);
}

void query(ll node,ll b,ll e,ll p){
    if(b==e){
        if(b==p)
            cur=segtree[node].val+segtree[node].x;
        return;
    }
    if(b>p||e<p) return;

    ll mid=(b+e)/2;
    ll left=2*node;
    ll right=left+1;

    if(segtree[node].f){
        segtree[left].f=true;
        segtree[left].x+=segtree[node].x;
        segtree[right].f=true;
        segtree[right].x+=segtree[node].x;

        segtree[node].x=0;
        segtree[node].f=false;
    }

    
    query(left,b,mid,p);
    query(right,mid+1,e,p);
    return;
}
*/
//Hashing
void  Hashing::build(){
    ull p=1,b=97;
    forwdhash[0]=0;
    for(ll i=1;i<=n;i++){
        forwdhash[i]=forwdhash[i-1]+a[i]*p;
        po[i]=p;
        p*=b;
    }

    p=1;
    revhash[n+1]=0;
    for(ll i=n;i>=1;i--){
        revhash[i]=revhash[i+1]+a[i]*p;
        revpo[i]=p;
        p*=b;
    }

}

bool Hashing::isPlain(ll l,ll r){
    ull x=forwdhash[r]-forwdhash[l-1];
    ull y=revhash[l]-revhash[r+1];
    ull rpo=n-r+1;
    ull lpo=l;
    if(lpo==rpo){
        if(x==y) return true;
        else return false;
    }
    else if(lpo>rpo){
        ull d=lpo-rpo;
        if(x==y*po[d+1]) return true;
        else return false;
    }
    else{
        ull d=rpo-lpo;
        if(x*po[d+1]==y) return true;
        else return false;
    }

}

ll Miscellaneous::lis(vector<ll> const& a) {
    //nlog(n)
    ll n = a.size();
    const ll INF = 1e9;
    vector<ll> d(n+1, INF);
    d[0] = -INF;
    for (ll i = 0; i < n; i++) {
        ll l = upper_bound(d.begin(), d.end(), a[i]) - d.begin();
        if (d[l-1] < a[i] && a[i] < d[l])
            d[l] = a[i];
    }
    ll ans = 0;
    for (ll l = 0; l <= n; l++) {
        if (d[l] < INF)
            ans = l;
    }
    return ans;
}

bool Miscellaneous::issubsequence(string& s1, string& s2)
{   //finds if s1 is a subsec of s2 or not
    ll n = s1.length(), m = s2.length();
    ll i = 0, j = 0;
    while (i < n && j < m) {
        if (s1[i] == s2[j])
            i++;
        j++;
    }
    return i == n;
}
dob Geometry::SQR(dob x){return x*x;}
dob Geometry::DOT(dob x1,dob y1,dob x2,dob y2){ return x1*x2+y1*y2;}
dob Geometry::CROSS(dob x1,dob y1,dob x2,dob y2){ return x1*y2-x2*y1;}
dob Geometry::PolarR(dob x,dob y){ return sqrt(SQR(x)+SQR(y));}
dob Geometry::PolarAngle(dob x,dob y){ return atan2(dob(y),dob(x));}
dob Geometry::CartesianX(dob r,dob ang){ return r*cos(dob(ang));}
dob Geometry::CartesianY(dob r,dob ang){ return r*sin(dob(ang));}
dob Geometry::CartesianDis(dob x1,dob y1,dob x2,dob y2){return sqrt( SQR(x1-x2) + SQR(y1-y2) );}
dob Geometry::PolarDis(dob r1,dob ang1,dob r2,dob ang2){return sqrt(SQR(r1)+SQR(r2)-2.00*r1*r2*cos(ang1-ang2));}
dob Geometry::MiddlePoint(dob x1,dob x2){return (x1+x2)/2;}
dob Geometry::InternalSection(dob x1,dob x2,dob m1,dob m2){ return (m1*x2+m2*x1)/(m1+m2);}
dob Geometry::ExternalSection(dob x1,dob x2,dob m1,dob m2){ return (m1*x2-m2*x1)/(m1-m2);}
dob Geometry::AreaOfTriangle(dob x1,dob y1,dob x2,dob y2,dob x3,dob y3){ return (x1*y2+x2*y3+x3*y1)-(y1*x2+y2*x3+y3*x1);}
dob Geometry::GradAngle(dob x1,dob y1,dob x2,dob y2){ return atan2( dob(y1-y2) , dob(x1-x2));}
dob Geometry::Grad(dob x1,dob y1,dob x2,dob y2){ return (x1-x2)==0.00?1.00*1e18:( dob(y1-y2) / dob(x1-x2));}
dob Geometry::Line(dob m,dob c,dob x,dob y){ return y-m*x-c;}
dob Geometry::Line(dob m,dob x1,dob y1,dob x,dob y){ return (y-y1)-m*(x-x1);}
dob Geometry::Line(dob x1,dob y1,dob x2,dob y2,dob x,dob y){ return (y2-y1)*x-(x2-x1)*y+(x1*y2-y1*x2);}
dob Geometry::OrthogonalLine(dob x1,dob y1,dob x2,dob y2,dob lx,dob ly,dob ox,dob oy){ return -lx*(x2-x1)-ly*(y2-y1)+(y2-y1)*ox+(x2-x1)*oy;}
dob Geometry::LineOrthogonalLineIntersecX(dob x1,dob y1,dob x2,dob y2,dob ox,dob oy){ return ((x1-x2)*((y2-y1)*ox+(x2-x1)*oy)-(y1-y2)*(x1*y2-x2*y1))/((y2-y1)*(y1-y2)-SQR(x1-x2));}
dob Geometry::LineOrthogonalLineIntersectY(dob x1,dob y1,dob x2,dob y2,dob ox,dob oy){ return ((x1-x2)*(x1*y2-x2*y1)-(y2-y1)*((y2-y1)*ox+(x2-x1)*oy))/((y2-y1)*(y1-y2)-SQR(x1-x2));}
bool Geometry::IsParallel(dob x1,dob y1,dob x2,dob y2){ (CROSS(x1,y1,x2,y2)==0.00)?true:false;}
bool Geometry::IsOrthogonal(dob x1,dob y1,dob x2,dob y2){ (DOT(x1,y1,x2,y2)==0.00)?true:false;}
dob Geometry::IntersecX(dob a1,dob b1,dob c1,dob a2,dob b2,dob c2){ return (b1*c2-b2*c1)/(a1*b2-a2*b1);}
dob Geometry::IntersecY(dob a1,dob b1,dob c1,dob a2,dob b2,dob c2){ return (b1*c2-b2*c1)/(a1*b2-a2*b1);}
dob Geometry::PointToLineDis(dob a,dob b,dob c,dob x,dob y){ return (a*x+b*y+c)/sqrt(SQR(a)+SQR(b));}
dob Geometry::PointToLineDis(dob x1,dob y1,dob x2,dob y2,dob x,dob y){ return PointToLineDis(y2-y1,x1-x2,x1*y2-x2*y1,x,y);}
dob Geometry::LineToLineDis(dob a,dob b,dob c1,dob c2){ return abs(c1-c2)/sqrt(SQR(a)+SQR(b));}
ll Geometry::isLeft(CartesianPoint a,CartesianPoint b,CartesianPoint c) { return (a.x-b.x)*(b.y-c.y)-(b.x-c.x)*(a.y-b.y);}

/*
The weapons in my hand:
**STL
*Number theory
**Geometry

*sqrt decomposition
*mo's algorithm
**segment tree
**Disjoint set union (DSU)
*trie/ prefix tree
*merge sort tree
*tree dp
*pbds to delete elements in logn time

**bfs
**dfs
*2d grid a bfs
*2d grid a dfs
*lca
**heavy light decomposition/sack/dsu on tree

*classical dp
*bitmask dp
*digit dp

*hashing
*matrix exponenciation
*cumulative sum
*sliding window technique.
*two pointers
*Greedy
*sparse table
*Binary Search & Binary search on answer

*/
