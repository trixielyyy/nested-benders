// Nested (node-based) Benders prototype for the three-stage stochastic multiple knapsack.
//
// Structure:
//   root master : stage-1 vars x,w (integer) + theta_lambda
//   node master : stage-2 vars y^lambda (integer) + psi_s per s in S(lambda)
//   subproblem  : stage-3 vars z^s (continuous LP)
//
// Two modes:
//
//   exact (default) : the node is a MIP, solved to optimality for the current x-hat with
//                     stage-3 cuts added between node solves. The cut sent to the root is
//                     an integer optimality cut: theta_lambda >= V_lambda(x-hat), active only
//                     when x = x-hat. The root therefore needs no convexity assumption.
//
//   --relaxed       : stage-2 integrality is dropped, the node is an LP, and the cut sent to
//                     the root is the usual subgradient cut. The bound converges to the LP
//                     relaxation, not to the MIP optimum.
//
// Usage: ./nested_benders <instance-prefix> [max_outer] [--relaxed]

#include "GurobiSolve3Stage.h"   // also pulls in ProblemParams via ProblemParams.cpp

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

double objOf(const ProblemParams& p, const std::string& scen, const std::string& var) {
    auto s = p.stoc_obj.find(scen);
    if (s == p.stoc_obj.end()) return 0.0;
    auto v = s->second.find(var);
    if (v == s->second.end()) return 0.0;
    return v->second * p.probs.at(scen);
}

double ubOf(const ProblemParams& p, const std::string& var) {
    return p.ub.count(var) ? p.ub.at(var) : GRB_INFINITY;
}

char senseOf(const ProblemParams& p, const std::string& row) {
    auto it = p.sense.find(row);
    if (it == p.sense.end()) return 'E';
    return it->second.empty() ? 'E' : it->second[0];
}

char gurobiSense(char s) {
    if (s == 'L') return GRB_LESS_EQUAL;
    if (s == 'G') return GRB_GREATER_EQUAL;
    return GRB_EQUAL;
}

double rhsOf(const ProblemParams& p, const std::string& row) {
    return p.rhs.count(row) ? p.rhs.at(row) : 0.0;
}

// The MPS marker alone does not prove a binary. A general integer keeps GRB_INTEGER,
// so the exact-mode guard below rejects it instead of building wrong link rows.
char varType(const ProblemParams& p, const std::string& var) {
    if (!(p.type.count(var) && p.type.at(var) == "I")) return GRB_CONTINUOUS;
    return ubOf(p, var) <= 1.0 + 1e-9 ? GRB_BINARY : GRB_INTEGER;
}

struct ExtensiveResult {
    double obj = std::nan("");
    int vars = 0;
    int rows = 0;
};

// ---------------------------------------------------------------------------
// Reference model: the full extensive form in one Gurobi model.
// ---------------------------------------------------------------------------
ExtensiveResult solveExtensive(GRBEnv& env, const ProblemParams& p, bool relaxStage2, double tl) {
    ExtensiveResult out;
    GRBModel m(env);
    m.set(GRB_IntParam_OutputFlag, 0);
    if (tl > 0) m.set(GRB_DoubleParam_TimeLimit, tl);

    std::map<std::string, GRBVar> x1;
    std::map<std::string, std::map<std::string, GRBVar>> y2, z3;
    for (const auto& kv : p.col_stage) {
        const std::string& n = kv.first;
        const std::string& st = kv.second;
        if (st == "STAGE-1") {
            x1[n] = m.addVar(0.0, ubOf(p, n), p.obj.count(n) ? p.obj.at(n) : 0.0,
                             varType(p, n), n);
        } else if (st == "STAGE-2") {
            for (const auto& s : p.stage2_scenarios)
                y2[s][n] = m.addVar(0.0, ubOf(p, n), objOf(p, s, n),
                                    relaxStage2 ? GRB_CONTINUOUS : varType(p, n), n + "_" + s);
        } else if (st == "STAGE-3") {
            for (const auto& s : p.stage3_scenarios)
                z3[s][n] = m.addVar(0.0, ubOf(p, n), objOf(p, s, n), GRB_CONTINUOUS, n + "_" + s);
        }
    }
    m.update();

    std::map<std::string, GRBLinExpr> r1;
    std::map<std::string, std::map<std::string, GRBLinExpr>> r2, r3;
    for (const auto& kv : p.row_stage) {
        if (kv.second == "STAGE-1") r1[kv.first] = GRBLinExpr();
        else if (kv.second == "STAGE-2") for (const auto& s : p.stage2_scenarios) r2[s][kv.first] = GRBLinExpr();
        else if (kv.second == "STAGE-3") for (const auto& s : p.stage3_scenarios) r3[s][kv.first] = GRBLinExpr();
    }
    for (const auto& vk : p.coeff) {
        const std::string& v = vk.first;
        const std::string& vs = p.col_stage.at(v);
        for (const auto& rk : vk.second) {
            const std::string& row = rk.first;
            double c = rk.second;
            const std::string& rs = p.row_stage.at(row);
            if (rs == "STAGE-1") {
                r1[row] += c * x1[v];
            } else if (rs == "STAGE-2") {
                for (const auto& s : p.stage2_scenarios) {
                    if (vs == "STAGE-1") r2[s][row] += c * x1[v];
                    else if (vs == "STAGE-2") r2[s][row] += c * y2[s][v];
                }
            } else if (rs == "STAGE-3") {
                for (const auto& s : p.stage3_scenarios) {
                    if (vs == "STAGE-1") {
                        r3[s][row] += c * x1[v];
                    } else if (vs == "STAGE-2") {
                        std::string par = p.parentStage2For(s);
                        if (!par.empty()) r3[s][row] += c * y2[par][v];
                    } else if (vs == "STAGE-3") {
                        r3[s][row] += c * z3[s][v];
                    }
                }
            }
        }
    }
    for (auto& kv : r1) m.addConstr(kv.second, gurobiSense(senseOf(p, kv.first)), rhsOf(p, kv.first), kv.first);
    for (auto& sk : r2) for (auto& kv : sk.second)
        m.addConstr(kv.second, gurobiSense(senseOf(p, kv.first)), rhsOf(p, kv.first), kv.first);
    for (auto& sk : r3) for (auto& kv : sk.second)
        m.addConstr(kv.second, gurobiSense(senseOf(p, kv.first)), rhsOf(p, kv.first), kv.first);
    m.optimize();
    if (m.get(GRB_IntAttr_Status) == GRB_OPTIMAL || m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
        out.obj = m.get(GRB_DoubleAttr_ObjVal);
        out.vars = m.get(GRB_IntAttr_NumVars);
        out.rows = m.get(GRB_IntAttr_NumConstrs);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Nested Benders
// ---------------------------------------------------------------------------
struct Stage3Cut {
    std::string s3;                       // stage-3 scenario
    double cst = 0.0;                     // psi_s >= cst + sum a_j x_j + sum b_j y_j
    std::map<std::string, double> a;      // coefficients on stage-1 variables
    std::map<std::string, double> b;      // coefficients on stage-2 variables
    GRBConstr con{};
};

struct Node {
    std::unique_ptr<GRBModel> model;
    std::map<std::string, GRBVar> y;
    std::map<std::string, GRBVar> psi;
    std::map<std::string, GRBConstr> rowCon;   // stage-2 rows, RHS follows x-hat
    std::vector<Stage3Cut> cuts;
    int cutsAdded = 0;
    int solves = 0;
    int innerCapHits = 0;
};

struct Stage3Sub {
    std::unique_ptr<GRBModel> model;
    std::map<std::string, GRBConstr> rowCon;
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <instance-prefix> [max_outer] [--relaxed]\n", argv[0]);
        return 2;
    }
    std::vector<std::string> args(argv + 1, argv + argc);
    bool relaxed = false;
    for (const auto& a : args) if (a == "--relaxed") relaxed = true;

    const std::string prefix = args[0];
    int maxOuter = 200;
    if (args.size() > 1 && args[1][0] != '-') maxOuter = std::atoi(args[1].c_str());
    const double tol = 1e-6;
    const double tl = 900.0;
    const int maxInner = 40;

    ProblemParams p;
    p.parse_tim(prefix + "_3stage.tim");
    p.parse_cor(prefix + "_3stage.cor");
    p.parse_sto(prefix + "_extended.sto");

    // The integer optimality cut needs binary stage-1 variables.
    if (!relaxed) {
        for (const auto& kv : p.col_stage)
            if (kv.second == "STAGE-1" && varType(p, kv.first) != GRB_BINARY) {
                std::printf("stage-1 variable %s is continuous, falling back to --relaxed\n", kv.first.c_str());
                relaxed = true;
                break;
            }
    }

    const auto& S2 = p.stage2_scenarios;
    const auto& S3 = p.stage3_scenarios;

    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);

    const auto t0 = std::chrono::steady_clock::now();
    ExtensiveResult extMip = solveExtensive(env, p, false, tl);
    ExtensiveResult extLp = solveExtensive(env, p, true, tl);
    const auto t1 = std::chrono::steady_clock::now();
    std::printf("mode        %s\n", relaxed ? "relaxed (LP nodes)" : "exact (integer nodes)");
    std::printf("extensive   MIP %.4f  vars %d  rows %d", extMip.obj, extMip.vars, extMip.rows);
    if (!std::isnan(extLp.obj)) std::printf("   LP %.4f", extLp.obj);
    std::printf("   time %.2f s\n", std::chrono::duration<double>(t1 - t0).count());

    // ---- root master ------------------------------------------------------
    GRBModel root(env);
    root.set(GRB_IntParam_OutputFlag, 0);
    root.set(GRB_DoubleParam_TimeLimit, tl);
    std::map<std::string, GRBVar> x1;
    for (const auto& kv : p.col_stage)
        if (kv.second == "STAGE-1")
            x1[kv.first] = root.addVar(0.0, ubOf(p, kv.first),
                                       p.obj.count(kv.first) ? p.obj.at(kv.first) : 0.0,
                                       varType(p, kv.first), kv.first);
    std::map<std::string, GRBVar> theta;
    for (const auto& l : S2)
        theta[l] = root.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "theta_" + l);
    // theta needs a finite lower bound. Without one the root gives every unvisited pattern -inf
    // and the loop never certifies. ponytail: 0 assumes nonnegative stage-2 and stage-3 costs.
    double thetaLb = 0.0;
    for (const auto& sk : p.stoc_obj)
        for (const auto& vk : sk.second)
            if (vk.second < 0.0) thetaLb = -GRB_INFINITY;
    for (auto& kv : theta) kv.second.set(GRB_DoubleAttr_LB, thetaLb);
    if (thetaLb < 0.0 && !relaxed) {
        std::printf("negative cost coefficients: the big-M value M = V is not valid, falling back to --relaxed\n");
        relaxed = true;
    }
    root.update();
    {
        std::map<std::string, GRBLinExpr> rows;
        for (const auto& kv : p.row_stage)
            if (kv.second == "STAGE-1") rows[kv.first] = GRBLinExpr();
        for (const auto& vk : p.coeff) {
            if (p.col_stage.at(vk.first) != "STAGE-1") continue;
            for (const auto& rk : vk.second)
                if (p.row_stage.at(rk.first) == "STAGE-1") rows[rk.first] += rk.second * x1[vk.first];
        }
        for (auto& kv : rows)
            root.addConstr(kv.second, gurobiSense(senseOf(p, kv.first)), rhsOf(p, kv.first), kv.first);
    }

    // ---- node masters -----------------------------------------------------
    const char yType = relaxed ? GRB_CONTINUOUS : GRB_BINARY;
    std::map<std::string, Node> nodes;
    for (const auto& l : S2) {
        Node& nd = nodes[l];
        nd.model = std::make_unique<GRBModel>(env);
        nd.model->set(GRB_IntParam_OutputFlag, 0);
        nd.model->set(GRB_DoubleParam_TimeLimit, tl);
        for (const auto& kv : p.col_stage)
            if (kv.second == "STAGE-2" && !kv.first.empty())
                nd.y[kv.first] = nd.model->addVar(0.0, ubOf(p, kv.first), objOf(p, l, kv.first),
                                                  yType, kv.first);
        std::map<std::string, GRBLinExpr> rows;
        for (const auto& kv : p.row_stage)
            if (kv.second == "STAGE-2") rows[kv.first] = GRBLinExpr();
        for (const auto& vk : p.coeff) {
            const std::string& v = vk.first;
            if (p.col_stage.at(v) != "STAGE-2") continue;
            for (const auto& rk : vk.second)
                if (p.row_stage.at(rk.first) == "STAGE-2") rows[rk.first] += rk.second * nd.y[v];
        }
        nd.model->update();
        for (auto& kv : rows)
            nd.rowCon[kv.first] = nd.model->addConstr(kv.second, gurobiSense(senseOf(p, kv.first)), 0.0, kv.first);
        for (const auto& s : S3)
            if (p.parentStage2For(s) == l)
                nd.psi[s] = nd.model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "psi_" + s);
        nd.model->update();
    }

    // ---- stage-3 subproblems ---------------------------------------------
    std::map<std::string, Stage3Sub> subs;
    for (const auto& s : S3) {
        Stage3Sub& sp = subs[s];
        sp.model = std::make_unique<GRBModel>(env);
        sp.model->set(GRB_IntParam_OutputFlag, 0);
        sp.model->set(GRB_DoubleParam_TimeLimit, tl);
        std::map<std::string, GRBVar> z;
        for (const auto& kv : p.col_stage)
            if (kv.second == "STAGE-3" && !kv.first.empty())
                z[kv.first] = sp.model->addVar(0.0, ubOf(p, kv.first), objOf(p, s, kv.first),
                                               GRB_CONTINUOUS, kv.first);
        std::map<std::string, GRBLinExpr> rows;
        for (const auto& kv : p.row_stage)
            if (kv.second == "STAGE-3") rows[kv.first] = GRBLinExpr();
        for (const auto& vk : p.coeff) {
            if (p.col_stage.at(vk.first) != "STAGE-3") continue;
            for (const auto& rk : vk.second)
                if (p.row_stage.at(rk.first) == "STAGE-3") rows[rk.first] += rk.second * z[vk.first];
        }
        sp.model->update();
        for (auto& kv : rows)
            sp.rowCon[kv.first] = sp.model->addConstr(kv.second, gurobiSense(senseOf(p, kv.first)),
                                                      rhsOf(p, kv.first), kv.first);
    }

    // ---- iteration --------------------------------------------------------
    std::map<std::string, double> xHat;
    double lb = -GRB_INFINITY, prevLb = -GRB_INFINITY;
    int outer = 0;
    long subSolves = 0;
    double lastNodeSum = 0.0;
    double bestFeasible = GRB_INFINITY;   // true cost of the last examined x-hat
    bool capHit = false;
    bool certified = false;

    for (; outer < maxOuter; ++outer) {
        root.optimize();
        if (root.get(GRB_IntAttr_Status) != GRB_OPTIMAL) { std::printf("root not optimal\n"); break; }
        lb = root.get(GRB_DoubleAttr_ObjVal);
        for (auto& kv : x1) xHat[kv.first] = kv.second.get(GRB_DoubleAttr_X);

        lastNodeSum = 0.0;
        for (const auto& l : S2) {
            Node& nd = nodes[l];

            for (auto& rk : nd.rowCon) {
                double rhs = rhsOf(p, rk.first);
                for (const auto& vk : p.coeff) {
                    if (p.col_stage.at(vk.first) != "STAGE-1") continue;
                    auto c = vk.second.find(rk.first);
                    if (c != vk.second.end()) rhs -= c->second * xHat[vk.first];
                }
                rk.second.set(GRB_DoubleAttr_RHS, rhs);
            }
            for (auto& cut : nd.cuts) {
                double rhs = cut.cst;
                for (const auto& a : cut.a) rhs += a.second * xHat[a.first];
                cut.con.set(GRB_DoubleAttr_RHS, rhs);
            }

            // Inner loop: solve the node master, then check the stage-3 cost at its solution.
            double V = 0.0;
            std::map<std::string, double> yHat;
            bool innerDone = false;
            for (int inner = 0; inner < maxInner; ++inner) {
                nd.model->optimize();
                ++nd.solves;
                const int st = nd.model->get(GRB_IntAttr_Status);
                if (st != GRB_OPTIMAL) {
                    std::printf("node %s status %d at inner iteration %d\n", l.c_str(), st, inner);
                    return 1;
                }
                V = nd.model->get(GRB_DoubleAttr_ObjVal);
                for (auto& kv : nd.y) yHat[kv.first] = kv.second.get(GRB_DoubleAttr_X);

                bool added = false;
                for (const auto& s : S3) {
                    if (p.parentStage2For(s) != l) continue;
                    Stage3Sub& sp = subs[s];
                    for (auto& rk : sp.rowCon) {
                        double rhs = rhsOf(p, rk.first);
                        for (const auto& vk : p.coeff) {
                            const std::string& vs = p.col_stage.at(vk.first);
                            auto c = vk.second.find(rk.first);
                            if (c == vk.second.end()) continue;
                            if (vs == "STAGE-1") rhs -= c->second * xHat[vk.first];
                            else if (vs == "STAGE-2" && p.parentStage2For(s) == l) rhs -= c->second * yHat[vk.first];
                        }
                        rk.second.set(GRB_DoubleAttr_RHS, rhs);
                    }
                    sp.model->optimize();
                    ++subSolves;
                    if (sp.model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                        std::printf("subproblem %s not optimal\n", s.c_str());
                        return 1;
                    }
                    const double phi = sp.model->get(GRB_DoubleAttr_ObjVal);
                    const double psiVal = nd.psi[s].get(GRB_DoubleAttr_X);
                    if (phi > psiVal + tol) {
                        Stage3Cut cut;
                        cut.s3 = s;
                        for (const auto& rk : sp.rowCon) {
                            const double pi = rk.second.get(GRB_DoubleAttr_Pi);
                            if (pi == 0.0) continue;
                            cut.cst += pi * rhsOf(p, rk.first);
                            for (const auto& vk : p.coeff) {
                                auto c = vk.second.find(rk.first);
                                if (c == vk.second.end()) continue;
                                const std::string& vs = p.col_stage.at(vk.first);
                                if (vs == "STAGE-1") cut.a[vk.first] -= pi * c->second;
                                else if (vs == "STAGE-2") cut.b[vk.first] -= pi * c->second;
                            }
                        }
                        GRBLinExpr e = GRBLinExpr(nd.psi[s]);
                        for (const auto& b : cut.b) e += -b.second * nd.y[b.first];
                        cut.con = nd.model->addConstr(e, GRB_GREATER_EQUAL, 0.0, "cut_" + s);
                        nd.cuts.push_back(cut);
                        ++nd.cutsAdded;
                        added = true;
                    }
                }
                if (!added) { innerDone = true; break; }
                nd.model->update();
            }
            if (!innerDone) { ++nd.innerCapHits; capHit = true; }
            lastNodeSum += V;

            if (relaxed) {
                // theta_lambda >= V + g . (x - x-hat), g = subgradient of the node LP value
                std::map<std::string, double> g;
                for (const auto& rk : nd.rowCon) {
                    const double pi = rk.second.get(GRB_DoubleAttr_Pi);
                    if (pi == 0.0) continue;
                    for (const auto& vk : p.coeff) {
                        if (p.col_stage.at(vk.first) != "STAGE-1") continue;
                        auto c = vk.second.find(rk.first);
                        if (c != vk.second.end()) g[vk.first] -= pi * c->second;
                    }
                }
                for (const auto& cut : nd.cuts) {
                    const double mu = cut.con.get(GRB_DoubleAttr_Pi);
                    if (mu == 0.0) continue;
                    for (const auto& a : cut.a) g[a.first] += mu * a.second;
                }
                GRBLinExpr e = GRBLinExpr(theta[l]);
                double rhs = V;
                for (const auto& kv : g) {
                    e += -kv.second * x1[kv.first];
                    rhs -= kv.second * xHat[kv.first];
                }
                root.addConstr(e, GRB_GREATER_EQUAL, rhs, "bcut_" + l);
            } else {
                // Integer optimality cut (Laporte and Louveaux 1993), big-M form:
                //   theta_lambda >= V - M * (number of stage-1 coordinates that differ from x-hat)
                // M = V is valid here: every cost is nonnegative, so V >= 0 and theta_lambda >= 0,
                // and one deviation already makes the row vacuous. No new binary is needed.
                GRBLinExpr dev;
                for (const auto& kv : xHat) {
                    const GRBVar& xv = x1[kv.first];
                    if (kv.second > 0.5) { dev += 1.0; dev -= xv; }
                    else                 { dev += xv; }
                }
                GRBLinExpr row = GRBLinExpr(theta[l]);
                row += V * dev;
                root.addConstr(row, GRB_GREATER_EQUAL, V, "icut_" + l + "_" + std::to_string(outer));
            }
            root.update();
        }

        double fx = 0.0;
        for (const auto& kv : x1)
            if (p.obj.count(kv.first)) fx += p.obj.at(kv.first) * xHat[kv.first];
        bestFeasible = fx + lastNodeSum;

        const double rel = std::fabs(bestFeasible) > 1.0 ? std::fabs(bestFeasible) : 1.0;
        // Exact mode: the root value is a lower bound on the optimum, bestFeasible is an upper
        // bound. Meeting in the middle certifies the incumbent.
        if (!relaxed && thetaLb > -GRB_INFINITY && lb >= bestFeasible - 1e-7 * rel) {
            certified = true;
            ++outer;
            break;
        }
        if (relaxed && outer > 0 && std::fabs(lb - prevLb) <= 1e-7 * rel) { ++outer; break; }
        prevLb = lb;
    }

    const double reported = relaxed ? lb : bestFeasible;

    const auto t2 = std::chrono::steady_clock::now();
    std::printf("nested      obj %.4f%s   outer %d   node solves %ld   subproblem solves %ld   time %.2f s\n",
                reported, certified ? " (certified)" : " (not certified)", outer, [&] { long n = 0; for (auto& kv : nodes) n += kv.second.solves; return n; }(),
                subSolves, std::chrono::duration<double>(t2 - t0).count());
    std::printf("masters     root vars %d rows %d   node vars %d rows %d (first node)\n",
                root.get(GRB_IntAttr_NumVars), root.get(GRB_IntAttr_NumConstrs),
                nodes.begin()->second.model->get(GRB_IntAttr_NumVars),
                nodes.begin()->second.model->get(GRB_IntAttr_NumConstrs));
    int cuts = 0;
    for (auto& kv : nodes) cuts += kv.second.cutsAdded;
    std::printf("cuts        stage-3 %d   inner cap hits %d%s\n", cuts, [&] { int n = 0; for (auto& kv : nodes) n += kv.second.innerCapHits; return n; }(),
                capHit ? "  (a node hit the inner cap, V is not confirmed)" : "");
    if (!std::isnan(extMip.obj))
        std::printf("check       vs extensive MIP %+.4f (%.4f%%)\n", reported - extMip.obj,
                    100.0 * (reported - extMip.obj) / std::fabs(extMip.obj));
    if (!std::isnan(extLp.obj))
        std::printf("            vs extensive LP  %+.4e\n", reported - extLp.obj);
    return 0;
}
