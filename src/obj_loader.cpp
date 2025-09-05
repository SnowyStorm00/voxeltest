#include "obj_loader.h"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <tuple>
#include <array>

using namespace DirectX;

namespace {
struct Key { int v; int vt; int vn; bool operator==(const Key& o) const { return v==o.v && vt==o.vt && vn==o.vn; } };
struct KeyHash { size_t operator()(const Key& k) const { return ((k.v*73856093) ^ (k.vt*19349663) ^ (k.vn*83492791)); } };
}

bool LoadOBJ(const std::wstring& path, ModelMesh& out){
    std::ifstream f(path);
    if(!f) return false;
    std::vector<XMFLOAT3> pos, nor;
    std::vector<std::array<Key,4>> faces;
    std::string line;
    while(std::getline(f, line)){
        if(line.size()<2) continue;
        std::istringstream ss(line);
        std::string tag; ss >> tag;
        if(tag == "v"){
            XMFLOAT3 p; ss >> p.x >> p.y >> p.z; pos.push_back(p);
        } else if(tag == "vn"){
            XMFLOAT3 n; ss >> n.x >> n.y >> n.z; nor.push_back(n);
        } else if(tag == "f"){
            std::vector<Key> poly;
            std::string tok;
            while(ss >> tok){
                int v=0, vt=0, vn=0;
                // parse formats v, v//vn, v/vt/vn, v/vt
                size_t s1 = tok.find('/');
                if(s1==std::string::npos){ v = std::stoi(tok); }
                else {
                    size_t s2 = tok.find('/', s1+1);
                    v = std::stoi(tok.substr(0,s1));
                    if(s2==std::string::npos){ vt = std::stoi(tok.substr(s1+1)); }
                    else {
                        if(s2 > s1+1){ vt = std::stoi(tok.substr(s1+1, s2-(s1+1))); }
                        if(s2+1 < tok.size()) vn = std::stoi(tok.substr(s2+1));
                    }
                }
                poly.push_back({v,vt,vn});
            }
            // triangulate fan
            for(size_t i=1;i+1<poly.size();++i){
                std::array<Key,4> tri{ poly[0], poly[i], poly[i+1], Key{0,0,0} };
                faces.push_back(tri);
            }
        }
    }
    // build indices/vertices with dedup
    out.vertices.clear(); out.indices.clear();
    std::unordered_map<Key, uint32_t, KeyHash> map;
    map.reserve(faces.size()*3);
    auto get = [&](const Key& k){
        auto it = map.find(k);
        if(it != map.end()) return it->second;
        XMFLOAT3 p{0,0,0}, n{0,1,0};
        int vi = (k.v>0? k.v-1 : (k.v<0? (int)pos.size()+k.v : 0));
        if(vi>=0 && vi<(int)pos.size()) p = pos[vi];
        int ni = (k.vn>0? k.vn-1 : (k.vn<0? (int)nor.size()+k.vn : -1));
        if(ni>=0 && ni<(int)nor.size()) n = nor[ni];
        uint32_t idx = (uint32_t)out.vertices.size();
        out.vertices.push_back({p,n});
        map.emplace(k, idx);
        return idx;
    };
    for(auto &tri : faces){
        out.indices.push_back(get(tri[0]));
        out.indices.push_back(get(tri[1]));
        out.indices.push_back(get(tri[2]));
    }
    return !out.indices.empty();
}
