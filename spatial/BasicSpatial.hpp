#pragma once

#include "SpatialBase.h"
#include <vector>
#include <set>
#include <list>

const int DIV = 10;
const int MAX = 1000;

namespace utec {
namespace spatial {

  
/**
 * BasicSpatial implementation - Matrix/Hash Division Approach
 * 
 * Almacena los Points en una matriz de dimensiones MAX/DIV x MAX/DIV, agrupando los puntos que 
 * caigan en el rango DIV a dicha posición en la matriz mediante hash_func (int coord). Para la búsqueda, se aplica
 * esta misma hash_func para obtener la posicion inicial en la matriz, realizando una busqueda lineal sobre su contenido
 * para obtener el punto con distancia mínima. En caso no se haya encontrado, se expande el área de busqueda, ejecutando 
 * una busqueda lineal sobre las posiciones adyacentes hasta encontrar la distancia mínima, retornando el punto encontrado.
 * A pesar de esto, se estima que los casos promedio resultarían igual en una complejidad de O(n). Mediante pruebas experimentales, 
 * se realizó un benchmark del running time para la búsqueda de 10000 queries aleatorias con 10000 elementos en la estructura, variando el tamaño de 
 * DIV entre pruebas. En estas pruebas, se observó que el running time formaba una parábola, los valores de división cercanos a 1 y a 1000 
 * tenian un running time bastante alto, mientras que los intermedios corrian considerablemente mas rápido, teniendo su menor runtime 
 * con un DIV de alrededor de 10. También se hicieron pruebas con menores tamaños de elementos en la estructura, con lo cual se observó que 
 * con una menor densidad de elementos se obtenia un mejor desempeño con un valor de división mayor; con 100 elementos en la estructura,
 * un valor de DIV alrededor de 100 obtenía el mejor runtime. Para que esta estructura desempeñe correctamente, se debe elegir cuidadosamente 
 * este valor de división. Conforme aumente la densidad de los datos, reducir el tamaño de DIV debería mejorar el desempeño de la búsqueda en esos casos. 
 * Otra optimización podría darse con el uso de matrices esparsas, con el fin de aliviar un poco el uso de memoria.
 */

template <typename Point>
class BasicSpatial : public SpatialBase<Point> {
  private:
    std::vector<std::vector<std::list<Point>>> GRID = 
      std::vector<std::vector<std::list<Point>>>((MAX/DIV)+1, std::vector<std::list<Point>>((MAX/DIV)+1, std::list<Point>()));
  
  public:
  BasicSpatial() {};

  int hash_func (int coord) {
    return (float(coord) / DIV);
  }

  void insert(const Point& new_point) override{
    int x = hash_func (new_point.get(0));
    int y = hash_func (new_point.get(1));

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > hash_func(MAX)) x = hash_func(MAX);
    if (y > hash_func(MAX)) y = hash_func(MAX);

    auto& a = GRID[x][y];  
    a.push_back(new_point);
  }

  void find_nearest(int maxx, int maxy, int minx, int miny, const Point& reference_, Point& nearest_, double& dist_, bool& fstset_){
    if (minx < 0) minx = 0;
    if (miny < 0) miny = 0;
    if (maxx > hash_func(MAX)) maxx = hash_func(MAX);
    if (maxy > hash_func(MAX)) maxy = hash_func(MAX);
    
    // nearest_ = Point({-1,-1});
    // fstset_ = true;
    // horizontal
    for (int curx = maxx; curx >= minx; curx--){
      if (curx < 0) curx = 0;
      if (maxy < 0) maxy = 0;
      if (curx > hash_func(MAX)) curx = hash_func(MAX);
      if (maxy > hash_func(MAX)) maxy = hash_func(MAX);
      if (miny < 0) miny = 0;
      if (miny > hash_func(MAX)) miny = hash_func(MAX);

      // top horizontal
      for (auto i : GRID[curx][maxy]){
        if (!fstset_){
          nearest_ = i;
          fstset_ = true;
          dist_ = i.distance(reference_);
        }
        else{
          if (i.distance(reference_) < dist_){
            nearest_ = i;
            dist_ = i.distance(reference_);
          }
        }
      }
      // bottom horizontal
      for (auto i : GRID[curx][miny]){
        if (!fstset_){
          nearest_ = i;
          fstset_ = true;
          dist_ = i.distance(reference_);
        }
        else{
          if (i.distance(reference_) < dist_){
            nearest_ = i;
            dist_ = i.distance(reference_);
          }
        }
      }
    }

    // // vertical
    // for (int cury = maxy-1; cury >= miny+1; cury--){
    //   // left vertical
    //   for (auto i : GRID[maxx][cury]){
    //     if (!fstset_){
    //       nearest_ = i;
    //       fstset_ = true;
    //       dist_ = i.distance(reference_);
    //     }
    //     else{
    //       if (i.distance(reference_) < dist_){
    //         nearest_ = i;
    //         dist_ = i.distance(reference_);
    //       }
    //     }
    //   }
    //   // right vertical
    //   for (auto i : GRID[minx][cury]){
    //     if (!fstset_){
    //       nearest_ = i;
    //       fstset_ = true;
    //       dist_ = i.distance(reference_);
    //     }
    //     else{
    //       if (i.distance(reference_) < dist_){
    //         nearest_ = i;
    //         dist_ = i.distance(reference_);
    //       }
    //     }
    //   }
    // }
  }

  // El punto de referencia no necesariamente es parte del dataset
  Point nearest_neighbor(const Point& reference) override{ 
    int x = hash_func (reference.get(0));
    int y = hash_func (reference.get(1));
    
    int curmaxx = x; int curmaxy = y; 
    int curminx = x; int curminy = y;

    Point nearest = Point({0, 0}); 
    double dist = __DBL_MAX__; 
    bool fstset = 0;

    while (!fstset){
      if (curmaxy > hash_func(MAX) && curmaxx > hash_func(MAX) && curminx < 0 && curminy < 0)
        return Point({-1, -1});

      find_nearest(curmaxx,curmaxy,curminx,curminy,reference,nearest,dist,fstset);
      curmaxx = (curmaxx+1);
      curmaxy = (curmaxy+1);
      curminx = (curminx-1);
      curminy = (curminy-1);

      if (fstset)
        find_nearest(curmaxx,curmaxy,curminx,curminy,reference,nearest,dist,fstset);
    }
    return nearest; 
  }
};

}  // namespace spatial
}  // namespace utec
