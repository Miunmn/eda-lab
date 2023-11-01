#include "SStree.h"
using namespace std;
const int MAX = 1e9 + 7;

SsInnerNode::SsInnerNode(size_t d)
{
    centroid = Point(d);
    radius = 0;
}
SsLeaf::SsLeaf(size_t d)
{
    centroid = Point(d);
    radius = 0;
}

NType SsNode::varianceAlongDirection(const std::vector<Point> &centroids, size_t direction) const
{
    NType media = 0;
    NType varianza = 0;
    for (const Point &centroid : centroids)
    {
        media = media + centroid[direction];
    }
    media = media / centroids.size();
    for (const Point &centroid : centroids)
    {
        varianza = varianza + (centroid[direction] - media) * (centroid[direction] - media);
    }
    varianza = varianza / centroids.size();
    return varianza;
}

size_t SsNode::minVarianceSplit(size_t coordinateIndex)
{

    size_t ans = Settings::m;
    NType min = MAX;
    vector<Point> points = getEntriesCentroids();
    for (size_t i = Settings::m; i < Settings::M - Settings::m; i++)
    {
        vector<Point> l(points.begin(), points.begin() + i + 1);
        vector<Point> r(points.begin() + i + 1, points.end());
        if (varianceAlongDirection(r, coordinateIndex) + varianceAlongDirection(l, coordinateIndex) <= min)
        {
            min = varianceAlongDirection(r, coordinateIndex) + varianceAlongDirection(l, coordinateIndex);
            ans = i;
        }
    }

    return ans;
}

std::vector<std::string> SsTree::kNNQuery(const Point &center, size_t k) const
{
    // Estructura de datos para almacenar los k vecinos más cercanos
    std::priority_queue<Pair, std::vector<Pair>, Comparator> nearestNeighbors;

    // Inicializa la cola de prioridad con los k primeros nodos (recorrido BFS)
    std::queue<SsNode *> bfsQueue;
    bfsQueue.push(root);

    while (!bfsQueue.empty() && nearestNeighbors.size() < k)
    {
        SsNode *currentNode = bfsQueue.front();
        bfsQueue.pop();

        NType distance = ::distance(currentNode->centroid, center);
        nearestNeighbors.push(Pair(currentNode->centroid, distance));

        // Verifica si el nodo tiene hijos y agrégalos a la cola de prioridad
        if (!currentNode->isLeaf())
        {
            SsInnerNode *innerNode = dynamic_cast<SsInnerNode *>(currentNode);
            for (SsNode *child : innerNode->children)
            {
                NType childDistance = ::distance(child->centroid, center);
                bfsQueue.push(child);

                // Asegúrate de no agregar nodos fuera del rango de los k vecinos más cercanos
                if (nearestNeighbors.size() < k || childDistance < nearestNeighbors.top().distance)
                {
                    nearestNeighbors.push(Pair(child->centroid, childDistance));
                    // Si superamos k vecinos, elimina el más lejano
                    if (nearestNeighbors.size() > k)
                    {
                        nearestNeighbors.pop();
                    }
                }
            }
        }
    }

    // Recopila los k vecinos más cercanos en una lista
    std::vector<std::string> result;
    while (!nearestNeighbors.empty())
    {
        result.push_back(nearestNeighbors.top().point.path);
        nearestNeighbors.pop();
    }

    // Invierte la lista para que los vecinos más cercanos estén en el orden correcto
    std::reverse(result.begin(), result.end());

    return result;
}

size_t SsNode::directionOfMaxVariance() const
{
    vector<Point> centroids = getEntriesCentroids();
    size_t ans = 0;
    NType maxVariance = numeric_limits<NType>::min();
    for (size_t i = 0; i < centroids[0].dim(); ++i)
    {
        if (varianceAlongDirection(centroids, i) > maxVariance)
        {
            maxVariance = varianceAlongDirection(centroids, i);
            ans = i;
        }
    }

    return ans;
}

size_t SsNode::findSplitIndex()
{
    size_t coordinateIndex = directionOfMaxVariance();
    sortEntriesByCoordinate(coordinateIndex);
    return minVarianceSplit(coordinateIndex);
}
vector<Point> SsInnerNode::getEntriesCentroids() const
{
    vector<Point> centroids;
    for (auto child : children)
    {
        centroids.push_back(child->centroid);
    }
    return centroids;
}

void SsInnerNode::sortEntriesByCoordinate(size_t coordinateIndex)
{
    sort(children.begin(), children.end(), [coordinateIndex](const SsNode *n1, const SsNode *n2)
         { return n1->centroid[coordinateIndex] < n2->centroid[coordinateIndex]; });
}

std::pair<SsNode *, SsNode *> SsInnerNode::split()
{
    vector<SsNode *> l(children.begin(), children.begin() + findSplitIndex() + 1);
    vector<SsNode *> r(children.begin() + findSplitIndex() + 1, children.end());

    SsInnerNode *ln = new SsInnerNode();
    SsInnerNode *rn = new SsInnerNode();
    ln->children = l;
    rn->children = r;

    for (SsNode *child : l)
    {
        child->parent = ln;
    }

    for (SsNode *child : r)
    {
        child->parent = rn;
    }
    ln->updateBoundingEnvelope();
    rn->updateBoundingEnvelope();

    return make_pair(ln, rn);
}

SsNode *SsInnerNode::findClosestChild(const Point &target) const
{
    NType m1 = MAX;
    SsNode *ans = nullptr;
    for (auto child : children)
    {
        if (::distance(child->centroid, target) < m1)
        {
            m1 = ::distance(child->centroid, target);
            ans = child;
        }
    }
    return ans;
}

void SsInnerNode::updateBoundingEnvelope()
{
    vector<Point> points = getEntriesCentroids();
    centroid = Point(points[0].dim());
    for (const Point &point : points)
    {
        for (size_t i = 0; i < centroid.dim(); ++i)
        {
            centroid[i] = centroid[i] + point[i];
        }
    }
    centroid = centroid / points.size();

    radius = 0;
    for (const Point &point : points)
    {
        NType distance = ::distance(centroid, point) + radius;
        if (distance > radius)
        {
            radius = distance;
        }
    }
}

vector<Point> SsLeaf::getEntriesCentroids() const
{
    return points;
}

void SsLeaf::sortEntriesByCoordinate(size_t coordinateIndex)
{
    sort(points.begin(), points.end(), [coordinateIndex](const Point &p1, const Point &p2)
         { return p1[coordinateIndex] < p2[coordinateIndex]; });
}

pair<SsNode *, SsNode *> SsLeaf::split()
{
    vector<Point> l(points.begin(), points.begin() + findSplitIndex() + 1);
    vector<Point> r(points.begin() + findSplitIndex() + 1, points.end());
    SsLeaf *ln = new SsLeaf();
    SsLeaf *rn = new SsLeaf();
    ln->points = l;
    rn->points = r;
    ln->updateBoundingEnvelope();
    rn->updateBoundingEnvelope();

    return make_pair(ln, rn);
}

void SsLeaf::updateBoundingEnvelope()
{
    vector<Point> points = getEntriesCentroids();
    centroid = Point(points[0].dim());
    radius = 0;
    for (const Point &point : points)
    {
        for (size_t i = 0; i < centroid.dim(); ++i)
        {
            centroid[i] = centroid[i] + point[i];
        }
    }
    centroid = centroid / points.size();
    for (const Point &point : points)
    {
        if (::distance(centroid, point) >= radius)
        {
            radius = ::distance(centroid, point);
        }
    }
}

SsNode *SsInnerNode::insert(const Point &point)
{
    SsNode *closestChild = findClosestChild(point);
    closestChild->updateBoundingEnvelope();
    return closestChild->insert(point);
}

SsNode *SsLeaf::insert(const Point &point)
{
    points.push_back(point);
    updateBoundingEnvelope();

    if (points.size() > Settings::M)
    {
        pair<SsNode *, SsNode *> splitNodes = split();
        SsInnerNode *inner = new SsInnerNode();
        inner->children.push_back(splitNodes.first);
        inner->children.push_back(splitNodes.second);
        inner->updateBoundingEnvelope();
        splitNodes.first->parent = inner;
        splitNodes.second->parent = inner;

        if (this->parent)
        {
            SsInnerNode *parent = dynamic_cast<SsInnerNode *>(this->parent);
            for (size_t i = 0; i < parent->children.size(); ++i)
            {
                if (parent->children[i] == this)
                {
                    parent->children[i] = inner;
                    break;
                }
            }
            parent->updateBoundingEnvelope();
        }
        inner->parent = this->parent;
        if (this->parent)
        {
            this->parent->updateBoundingEnvelope();
        }
        return inner;
    }
    else
    {
        if (this->parent)
        {
            this->parent->updateBoundingEnvelope();
        }
        return this;
    }
}

void SsTree::insert(const Point &point)
{
    if (!root)
    {
        root = new SsLeaf();
        dynamic_cast<SsLeaf *>(root)->points.push_back(point);
        root->updateBoundingEnvelope();
        root->parent = nullptr;
    }
    else
    {
        SsNode *in = root->insert(point);
        if (in->parent == nullptr)
        {
            root = in;
        }
    }
}

void SsTree::insert(Point &point, const std::string &path)
{
    point.path = path;
    if (!root)
    {
        root = new SsLeaf();
        dynamic_cast<SsLeaf *>(root)->points.push_back(point);
        root->updateBoundingEnvelope();
        root->parent = nullptr;
    }
    else
    {
        SsNode *in = root->insert(point);
        if (in->parent == nullptr)
        {
            root = in;
        }
    }
}

void SsTree::build(const std::vector<Point> &points)
{
    for (const Point &point : points)
    {
        insert(point);
    }
}

bool SsNode::test(bool isRoot) const
{
    size_t count = 0;
    if (this->isLeaf())
    {
        const SsLeaf *leaf = dynamic_cast<const SsLeaf *>(this);
        count = leaf->points.size();

        for (const Point &point : leaf->points)
        {
            if (distance(this->centroid, point) > this->radius)
            {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    }
    else
    {
        const SsInnerNode *inner = dynamic_cast<const SsInnerNode *>(this);
        count = inner->children.size();

        for (const SsNode *child : inner->children)
        {
            if (distance(this->centroid, child->centroid) > this->radius)
            {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            if (!child->test(false))
            {
                return false;
            }
        }
    }

    if (!isRoot && (count < Settings::m || count > Settings::M))
    {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    if (!isRoot && !parent)
    {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}

void SsTree::test() const
{
    bool result = root->test(true);

    if (root->parent)
    {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result)
    {
        std::cout << "SS-Tree is valid!" << std::endl;
    }
    else
    {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}

void SsNode::print(size_t indent = 0) const
{
    for (size_t i = 0; i < indent; ++i)
    {
        std::cout << "  ";
    }

    // Imprime MAXormación del nodo.
    std::cout << "Centroid: " << centroid << ", Radius: " << radius;
    if (isLeaf())
    {
        const SsLeaf *leaf = dynamic_cast<const SsLeaf *>(this);
        std::cout << ", Points: [ ";
        for (const Point &p : leaf->points)
        {
            std::cout << p << " ";
        }
        std::cout << "]";
    }
    else
    {
        std::cout << std::endl;
        const SsInnerNode *inner = dynamic_cast<const SsInnerNode *>(this);
        for (const SsNode *child : inner->children)
        {
            child->print(indent + 1);
        }
    }
    std::cout << std::endl;
}
void SsTree::print() const
{
    if (root)
    {
        root->print();
    }
    else
    {
        std::cout << "Empty tree." << std::endl;
    }
}

void SsLeaf::saveToStream(std::ostream &out) const
{
    // Guardar centroid
    centroid.saveToFile(out, Settings::D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char *>(&radius_), sizeof(radius_));

    // Guardar el numero de puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char *>(&numPoints), sizeof(numPoints));

    // Guardar los puntos
    for (const auto &point : points)
    {
        point.saveToFile(out, Settings::D);
    }

    // Guardar las rutas (paths)
    size_t numPaths = paths.size();
    out.write(reinterpret_cast<const char *>(&numPaths), sizeof(numPaths));
    for (const auto &p : paths)
    {
        size_t pathLength = p.size();
        out.write(reinterpret_cast<const char *>(&pathLength), sizeof(pathLength));
        out.write(p.c_str(), (long)pathLength);
    }
}

void SsInnerNode::saveToStream(std::ostream &out) const
{
    // Guardar centroid
    centroid.saveToFile(out, Settings::D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char *>(&radius_), sizeof(radius_));

    // Guardar si apunta a nodos hoja
    bool pointsToLeafs = children[0]->isLeaf();
    out.write(reinterpret_cast<const char *>(&pointsToLeafs), sizeof(pointsToLeafs));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char *>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto &child : children)
    {
        child->saveToStream(out);
    }
}

void SsInnerNode::loadFromStream(std::istream &in)
{
    // Leer centroid
    centroid.readFromFile(in, Settings::D);

    // leer el valor del radio
    float radius_ = 0;
    in.read(reinterpret_cast<char *>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // leer si apunta a hojas o nodos internos
    bool pointsToLeaf = false;
    in.read(reinterpret_cast<char *>(&pointsToLeaf), sizeof(pointsToLeaf));

    // leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char *>(&numChildren), sizeof(numChildren));

    // leer hijos
    for (size_t i = 0; i < numChildren; ++i)
    {
        SsNode *child = pointsToLeaf ? static_cast<SsNode *>(new SsLeaf(Settings::D)) : static_cast<SsNode *>(new SsInnerNode(Settings::D));
        child->loadFromStream(in);
        children.push_back(child);
    }
}

void SsLeaf::loadFromStream(std::istream &in)
{
    // Leer centroid
    centroid.readFromFile(in, Settings::D);

    // Leer radio
    float radius_ = 0;
    in.read(reinterpret_cast<char *>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // Leer numero de puntos
    size_t numPoints;
    in.read(reinterpret_cast<char *>(&numPoints), sizeof(numPoints));

    // Leer puntos
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i)
    {
        points[i].readFromFile(in, Settings::D);
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char *>(&numPaths), sizeof(numPaths));
    paths.resize(numPaths);
    for (size_t i = 0; i < numPaths; ++i)
    {
        size_t pathLength;
        in.read(reinterpret_cast<char *>(&pathLength), sizeof(pathLength));
        char *buffer = new char[pathLength + 1];
        in.read(buffer, (long)pathLength);
        buffer[pathLength] = '\0';
        paths[i] = std::string(buffer);
        delete[] buffer;
    }
}

void SsTree::saveToFile(const std::string &filename) const
{
    std::ofstream out(filename, std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open file for writing");
    }

    // Guardar las dimensiones de la estructura
    out.write(reinterpret_cast<const char *>(&D), sizeof(D));

    // Guardar si el root es hija o nodo interno
    bool isLeaf = root->isLeaf();
    out.write(reinterpret_cast<const char *>(&isLeaf), sizeof(isLeaf));

    // Guardar el resto de la estructura
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in)
    {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root)
    {
        delete root;
        root = nullptr;
    }

    // Aquí se asume que el primer valor determina las dimensiones
    in.read(reinterpret_cast<char *>(&D), sizeof(D));

    // El segundo valor determina si el root es hoja
    bool isLeaf;
    in.read(reinterpret_cast<char *>(&isLeaf), sizeof(isLeaf));
    if (isLeaf)
    {
        root = new SsLeaf(D);
    }
    else
    {
        root = new SsInnerNode(D);
    }
    root->loadFromStream(in);
    in.close();
}
