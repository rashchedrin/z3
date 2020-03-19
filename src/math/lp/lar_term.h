/*++
  Copyright (c) 2017 Microsoft Corporation

  Module Name:

  <name>

  Abstract:

  <abstract>

  Author:

  Lev Nachmanson (levnach)

  Revision History:


  --*/
#pragma once
#include "math/lp/indexed_vector.h"
#include "util/map.h"

namespace lp {

class lar_term {
    // the term evaluates to sum of m_coeffs 
    typedef unsigned lpvar;
    map<tv, mpq, tv_hash, tv_equal > m_coeffs;
public:
    lar_term() {}
    void add_monomial(const mpq& c, tv j) {
        if (c.is_zero())
            return;
        auto *e = m_coeffs.find_core(j);        
        if (e == nullptr) {
            m_coeffs.insert(j, c);
        } else {
            e->get_data().m_value += c;
            if (e->get_data().m_value.is_zero())
                m_coeffs.erase(j);
        }
    }
    
    void add_var(tv j) {
        rational c(1);
        add_monomial(c, j);
    }

    bool is_empty() const {
        return m_coeffs.empty(); // && is_zero(m_v);
    }
    
    unsigned size() const { return static_cast<unsigned>(m_coeffs.size()); }

    template <typename T>
    const T & coeffs() const {
        return m_coeffs;
    }
    
    lar_term(const vector<std::pair<mpq, tv>>& coeffs) {
        for (const auto & p : coeffs) {
            add_monomial(p.first, p.second);
        }
    }
    bool operator==(const lar_term & a) const {  return false; } // take care not to create identical terms
    bool operator!=(const lar_term & a) const {  return ! (*this == a);}
    // some terms get used in add constraint
    // it is the same as the offset in the m_constraints

    vector<std::pair<mpq, tv>> coeffs_as_vector() const {
        vector<std::pair<mpq, tv>> ret;
        for (const auto & p :  m_coeffs) {
            ret.push_back(std::make_pair(p.m_value, p.m_key));
        }
        return ret;
    }

    // j is the basic variable to substitute
    void subst(tv j, indexed_vector<mpq> & li) {
        auto* it = m_coeffs.find_core(j);
        if (it == nullptr) return;
        const mpq & b = it->get_data().m_value;
        for (unsigned it_j :li.m_index) {
            add_monomial(- b * li.m_data[it_j], tv::var(it_j));
        }
        m_coeffs.erase(j);
    }

    // the monomial ax[j] is substituted by ax[k]
    void subst_index(tv j, tv k) {
        auto* it = m_coeffs.find_core(j);
        if (it == nullptr) return;
        mpq b = it->get_data().m_value;
        m_coeffs.erase(j);
        m_coeffs.insert(k, b);
    }
    
    bool contains(tv j) const {
        return m_coeffs.contains(j);
    }

    void negate() {
        for (auto & t : m_coeffs)
            t.m_value.neg();
    }

    template <typename T>
    T apply(const vector<T>& x) const {
        T ret(0);
        for (const auto & t : m_coeffs) {
            ret += t.m_value * x[t.m_key];
        }
        return ret;
    }
   
    void clear() {
        m_coeffs.reset();
    }

    class ival {
        tv m_tv;
        const mpq & m_coeff;
    public:
        ival(tv p, const mpq & val) : m_tv(p), m_coeff(val) { }
        tv tv() const { return m_tv; }
        const mpq & coeff() const { return m_coeff; }
    };
    
    class const_iterator {
        map<tv, mpq, tv_hash, tv_equal>::iterator m_it;
    public:
        ival operator*() const { return ival(m_it->m_key, m_it->m_value); }        
        const_iterator operator++() { const_iterator i = *this; m_it++; return i;  }
        const_iterator operator++(int) { m_it++; return *this; }
        const_iterator(map<tv, mpq, tv_hash, tv_equal>::iterator it) : m_it(it) {}
        bool operator==(const const_iterator &other) const { return m_it == other.m_it; }
        bool operator!=(const const_iterator &other) const { return !(*this == other); }
    };
   
    bool is_normalized() const {
        lpvar min_var = UINT_MAX;
        mpq c;
        for (const auto & p : *this) {
            SASSERT((p.tv().is_var()));
            if (p.tv().index() < min_var) {
                min_var = p.tv().index();
            }
        }
        lar_term r;
        for (const auto & p : *this) {
            if (p.tv().index() == min_var) {
                return p.coeff().is_one();
            }
        }
        lp_unreachable();
        return false;        
    }

    // a is the coefficient by which we divided the term to normalize it
    lar_term get_normalized_by_min_var(mpq& a) const {
        a = m_coeffs.begin()->m_value;
        if (a.is_one()) {
            return *this;
        }
        lar_term r;
        auto it = m_coeffs.begin();
        r.add_var(it->m_key);
        it++;
        for(;it != m_coeffs.end(); it++) {
            r.add_monomial(it->m_value / a, it->m_key);
        }
        return r;        
    }
    const_iterator begin() const { return m_coeffs.begin();}
    const_iterator end() const { return m_coeffs.end(); }
};
}
