#ifndef BNBTSPSOLVER_HPP
#define BNBTSPSOLVER_HPP

#include <iostream>
#include <cstdint>
#include <functional>
#include <vector>
#include <queue>
#include <stack>
#include <utility>

struct PartialSolution
{
public:
    friend class BnBTSPSol;
    PartialSolution(size_t n, const std::vector<std::vector<uint32_t> >& D) : n(n),Constraints(n*n,0)
    {
        Matrix=D;
        Reduced.resize(n*n);
        for(size_t i=0;i<n;++i)
        {
            for(size_t k=0;k<n;++k)
            {
                Reduced.at(i*n+k)=Matrix[i][k];
                if(i==k)
                {
                    Reduced.at(i*n+k)=INF;
                    Constraints.at(i*n+k)=-1;
                }
            }
        }
        Cost = 0;

        Reduce();
    }

    PartialSolution()
    {
    }

    bool operator>(const PartialSolution& other) const
    {
        return LowerBoundTimesTwo > other.LowerBoundTimesTwo;
    }
private:
    uint32_t INF=std::numeric_limits<uint32_t>::max();
    typedef std::pair < size_t, size_t > Edge;
    Edge NullEdge{-1,-1};
    enum class EdgeType
    {
        Outgoing,
        Incoming
    };

    PartialSolution WithEdge(Edge pivot)
    {
        auto i = pivot.first, j = pivot.second;

        PartialSolution child = *this;
        child.Cost += Matrix[i][j];
        for (size_t k = 0; k < n; k++)
        {
            child.Constraints[i*n+k] = child.Constraints[k*n+j] = -1;
            child.Reduced[i*n+k] = child.Reduced[k*n+j] = INF;
        }
        child.EnabledEdges++;
        child.Constraints[i*n+j] = 1;
        child.Constraints[j*n+i] = -1;

        auto subpathTo = child.TraverseSubPath(i, EdgeType::Outgoing);
        auto subpathFrom = child.TraverseSubPath(i, EdgeType::Incoming);
        if (subpathTo.size() + subpathFrom.size() - 1 != n)
        {
            size_t index=Index(subpathTo.back(), subpathFrom.back());
            child.Constraints[index] = -1;
            child.Reduced[index] = INF;
        }

        child.Reduce();
        return child;
    }

    PartialSolution WithoutEdge(Edge pivot)
    {
        auto i = pivot.first, j = pivot.second;

        PartialSolution child = *this;
        child.DisabledEdges++;
        child.Constraints[i*n+j] = -1;
        child.Reduced[i*n+j] = INF;
        child.Reduce(EdgeType::Outgoing, i);
        child.Reduce(EdgeType::Incoming, j);

        return child;
    }

    Edge ChoosePivotEdge(int)
    {
        auto minStride = [&](size_t except, size_t k, size_t kStride) {uint32_t m = INF; for (size_t i = 0; i < n; i++) if(i != except) m = std::min(m, Reduced[IK(i, k, kStride)]); return m;  };
        auto rowMin = [&](size_t k, size_t except) {return minStride(except, k, n); };
        auto columnMin = [&](size_t k, size_t except) {return minStride(except, k, 1); };

        uint32_t bestIncrease = 0;
        Edge bestPivot = NullEdge;
        for (size_t j = 0; j < n; j++)
        {
            if (Constraints[Index(0,j)] == 0 && Reduced[Index(0,j)] == 0)
            {
                auto increase = rowMin(0, j) + columnMin(j, 0);
                if (increase > bestIncrease)
                {
                    bestIncrease = increase;
                    bestPivot = Edge(0, j);
                }
            }
        }

        return bestPivot;
    }

    Edge ChoosePivotEdge()
    {
        auto minStride = [&](size_t except, size_t k, size_t kStride) {uint32_t m = INF; for (size_t i = 0; i < n; i++) if(i != except) m = std::min(m, Reduced[IK(i, k, kStride)]); return m;  };
        auto rowMin = [&](size_t k, size_t except) {return minStride(except, k, n); };
        auto columnMin = [&](size_t k, size_t except) {return minStride(except, k, 1); };

        uint32_t bestIncrease = 0;
        Edge bestPivot = NullEdge;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                if (Constraints[Index(i,j)] == 0 && Reduced[Index(i,j)] == 0)
                {
                    auto increase = rowMin(i, j) + columnMin(j, i);
                    if (increase > bestIncrease)
                    {
                        bestIncrease = increase;
                        bestPivot = Edge(i, j);
                    }
                }
            }
        }

        return bestPivot;
    }

    std::vector<size_t> TraverseSubPath(size_t cur, EdgeType edgeType)
    {
        auto stride = edgeType == EdgeType::Outgoing ? 1 : n;
        std::vector<size_t> subpath{ cur };
        for (size_t k = 0; k < n; k++)
        {
            size_t next = n;
            for (size_t i = 0; i < n; i++)
            {
                if (Constraints[IK(cur, i, stride)] == 1)
                {
                    next = i;
                    break;
                }
            }

            if (next == n)
                break;

            subpath.push_back(next);
            cur = next;
        }
        return subpath;
    }

    void Reduce(EdgeType edgeType, size_t i)
    {
        auto kStride = edgeType == EdgeType::Outgoing ? 1 : n;

        uint32_t min_element = INF;
        for (size_t k = 0; k < n; k++)
            if (Constraints[IK(i, k, kStride)] != -1)
                min_element = std::min(min_element, Reduced[IK(i, k, kStride)]);

        if (min_element < INF)
        {
            for (size_t k = 0; k < n; k++)
                Reduced[IK(i, k, kStride)] -= min_element;
            LowerBoundTimesTwo += min_element;
        }
    }

    void Reduce()
    {
        for (size_t i = 0; i < n; i++)
            Reduce(EdgeType::Outgoing, i);

        for (size_t j = 0; j < n; j++)
            Reduce(EdgeType::Incoming, j);
    }

    void Print()
    {
        std::cout << Cost - Matrix[n - 1][0] << std::endl;
        for (size_t i = 1; i < n; i++)
            std::cout << Path[i - 1] << " " << Path[i] << " " << Matrix[Path[i - 1]][Path[i]] << std::endl;
    }

    bool IsComplete()
    {
        Path = TraverseSubPath(0, EdgeType::Outgoing);
        //return Path.size() == n + 1 && Path[n - 1] == n - 1;
        return Path.size() == n + 1;
    }

    size_t IK(size_t i, size_t k, size_t kStride)
    {
        return (n + 1 - kStride)*i + kStride*k;
    }

    size_t Index(size_t row, size_t col)
    {
        return row*n+col;
    }

private:
    size_t EnabledEdges = 0, DisabledEdges = 0;
    size_t n = 0;
    uint32_t Cost = INF;
    uint32_t LowerBoundTimesTwo = 0;
    std::vector<std::vector<uint32_t> > Matrix;
    std::vector<uint32_t> Reduced;
    std::vector<int> Constraints;
    std::vector<size_t> Path;
};

class BnBTSPSol{
public:
    static void branch_and_bound(size_t n, std::vector<std::vector<uint32_t> > D)
    {
        PartialSolution bestCompleteSolution;
        //PartialSolution root = PartialSolution(n, D).WithEdge(Edge(n - 1, 0), D);
        auto rootPivot = PartialSolution(n, D).ChoosePivotEdge(int());
        PartialSolution root = PartialSolution(n, D).WithEdge(rootPivot);

        static std::priority_queue<PartialSolution, std::vector<PartialSolution>, std::greater<PartialSolution> > right;
        static std::stack<PartialSolution> left;
        left.push(root);

        while (!left.empty() || !right.empty())
        {
            auto currentSolution = !left.empty() ? left.top() : right.top();
            if (!left.empty())
                left.pop();
            else
                right.pop();

            if (currentSolution.IsComplete())
            {
                if (currentSolution.Cost < bestCompleteSolution.Cost)
                {
                    bestCompleteSolution = currentSolution;
                    //bestCompleteSolution.Print(D);
                }
            }
            else if (currentSolution.LowerBoundTimesTwo < bestCompleteSolution.Cost)
            {
                auto pivot = currentSolution.ChoosePivotEdge();
                if (pivot != currentSolution.NullEdge)
                {
                    auto withPivot = currentSolution.WithEdge(pivot);
                    auto withoutPivot = currentSolution.WithoutEdge(pivot);

                    if (withPivot.LowerBoundTimesTwo < bestCompleteSolution.Cost)
                        left.push(withPivot);

                    if (withoutPivot.LowerBoundTimesTwo < bestCompleteSolution.Cost)
                        right.push(withoutPivot);
                }
            }
        }
        bestCompleteSolution.Print();
    }
};


#endif // BNBTSPSOLVER_H
