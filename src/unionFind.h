// a class for transitive closure
class UnionFind{
public:
	UnionFind(int size);
	void makeSameGroup(int i,int j);

	int groupOf(int i);
	int nGroups();
	
	void returnGroups( std::vector< std::vector<int> > &groups );
	
private:
	std::vector<int> boss;
	void compressIfNeeded();
	std::vector<int> remap;
	int nGro;
	bool isCompressed ;
	int finalBoss(int i);
	
};

void UnionFind::returnGroups( std::vector< std::vector<int> > &groups){
	compressIfNeeded();
	groups.clear();
	groups.resize( nGroups() , std::vector<int>() ); // all groups start empty
	
	for (int i=0; i<(int)boss.size(); i++) {
		groups[ groupOf(i) ].push_back(i);
	}
}

UnionFind::UnionFind(int size){
	boss.resize(size);
	remap.resize(size);
	for (int i=0; i<(int)boss.size(); i++) {
		boss[i] = i;
		remap[i] = i;
	}
	nGro = size;
	isCompressed = true;
}

int UnionFind::nGroups(){
	compressIfNeeded();
	return nGro;
}

void UnionFind::compressIfNeeded(){
	if (isCompressed) return;
	nGro = 0;
	for (int i=0; i<(int)boss.size(); i++){
		if (boss[i]==i) remap[i] = nGro++;
		else remap[i] = -1;
	}
	isCompressed = true;
}

int UnionFind::finalBoss(int i){
	if (boss[i]==i) return i;
	boss[i] = boss[ boss[ i ] ]; // path shortening
	return ( finalBoss( boss[ i ] ) );
}

int UnionFind::groupOf(int i){
	compressIfNeeded();
	return remap[ finalBoss(i) ];
}

void UnionFind::makeSameGroup(int i, int j){
	int bi = finalBoss(i);
	int bj = finalBoss(j);
	if (bi == bj) return;
	if (rand()%2) boss[bi] = boss[bj]; else boss[bj] = boss[bi];
	isCompressed = false;
}