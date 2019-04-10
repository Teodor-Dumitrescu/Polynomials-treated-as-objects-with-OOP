#include <iostream>
#include <cstdlib>

using namespace std;

int prim(int x)
{
    int p=1;
    for(int i=2;i<=x/2;i++)
        if(x%i==0)
        p=0;
    return p;
}

class Monom
{
private:
    int grad;
    float coef;
public:
    int get_grad(){return grad;};
    void set_grad(int x){grad=x;};
    int get_coef(){return coef;};
    void set_coef(int x){coef=x;};
};

class Polinom
{
protected:
    int nr_monoame;
    Monom *m;
public:
    Polinom()
    {
        nr_monoame=1;
        m=new Monom;
        (*m).set_grad(0);
        (*m).set_coef(0);
    }
    Polinom(int nr)
    {
        nr_monoame=nr;
        m=new Monom[nr];
    }
    ~Polinom()
    {
        delete[] m;
    }
    Polinom(const Polinom& p)
    {
        nr_monoame=p.nr_monoame;
        m=new Monom[nr_monoame];
        for(int i=0;i<nr_monoame;i++)
        {
            m[i].set_grad(p.m[i].get_grad());
            m[i].set_coef(p.m[i].get_coef());
        }
    }
    Polinom& operator=(const Polinom& p)
    {
        if(this!=&p)
        {
            delete[] m;
            nr_monoame=p.nr_monoame;
            m=new Monom[nr_monoame];
            for(int i=0;i<nr_monoame;i++)
            {
                m[i].set_grad(p.m[i].get_grad());
                m[i].set_coef(p.m[i].get_coef());
            }
        }
        return *this;
    }
    void initializare();
    int crit_eisenstein();
    virtual void afisare(ostream &out)=0;
    friend istream& operator>>(istream& in, Polinom &p);
    friend ostream& operator<<(ostream& out, Polinom &p);
};

class Polinom_ireductibil: public Polinom
{
    public:
    Polinom_ireductibil(){};
    Polinom_ireductibil(int nr):Polinom(nr){};
    ~Polinom_ireductibil(){};
    Polinom_ireductibil(const Polinom_ireductibil& p):Polinom(p){};
    void afisare(ostream &out);
};

class Polinom_reductibil: public Polinom
{
    public:
    Polinom_reductibil(){};
    Polinom_reductibil(int nr):Polinom(nr){};
    ~Polinom_reductibil(){};
    Polinom_reductibil(const Polinom_reductibil& p):Polinom(p){};
    void afisare(ostream &out);
    Polinom_reductibil& horner(int x, int& bun);
};
Polinom_reductibil& Polinom_reductibil::horner(int x, int& bun)
{
    Polinom_reductibil *p=new Polinom_reductibil(m[0].get_grad()+1);
    int nr_monoame1=0;
    float vector1[m[0].get_grad()+1];
    for(int i=0;i<m[0].get_grad()+1;i++)
            vector1[i]=0;
    for(int i=0;i<nr_monoame;i++)
            vector1[m[0].get_grad()-m[i].get_grad()]=m[i].get_coef();
    float vector2[m[0].get_grad()+1];
    vector2[0]=vector1[0];
    for(int i=1;i<m[0].get_grad()+1;i++)
        vector2[i]=x*vector2[i-1]+vector1[i];
    if(vector2[m[0].get_grad()]==0)
        bun=1;
    //for(int i=0;i<m[0].get_grad()+1;i++)
        //cout<<i<<"   "<<vector2[i]<<endl;
    for(int i=0;i<m[0].get_grad()+1;i++)
    {
        if(vector2[i]!=0)
        {
            (*p).m[nr_monoame1].set_coef(vector2[i]);
            (*p).m[nr_monoame1].set_grad(m[0].get_grad()-i-1);
            nr_monoame1++;
        }
    }
    (*p).nr_monoame=nr_monoame1;
    return *p;
}
void Polinom_ireductibil::afisare(ostream &out)
{
    int i,grad;
    float coef;
    for(i=0;i<nr_monoame;i++)
    {
        grad=m[i].get_grad();
        coef=m[i].get_coef();
        if(i!=nr_monoame-1)
            out<<coef<<"x^"<<grad<<"+";
        else
            out<<coef<<"x^"<<grad;
    }
}
void Polinom_reductibil::afisare(ostream &out)
{
    int i,j,grad,bun=0;
    float coef;
    Polinom_reductibil p;
    for(i=0;i<nr_monoame;i++)
    {
        grad=m[i].get_grad();
        coef=m[i].get_coef();
        if(i!=nr_monoame-1)
            out<<coef<<"x^"<<grad<<"+";
        else
            out<<coef<<"x^"<<grad;
    }
    out<<endl;
    int x;
    if(m[nr_monoame-1].get_grad()>0)
    {
        p=(*this).horner(0,bun);
        out<<"x*(";
        for(j=0;j<nr_monoame;j++)
        {
            grad=p.m[j].get_grad();
            coef=p.m[j].get_coef();
            if(j!=p.nr_monoame-1)
                out<<coef<<"x^"<<grad<<"+";
            else
                out<<coef<<"x^"<<grad;
        }
        out<<")"<<endl;
    }
    else
    {
        for(x=-abs(m[nr_monoame-1].get_coef());x<=abs(m[nr_monoame-1].get_coef());x++)
        {
            if(x==0)
                continue;
            if(m[nr_monoame-1].get_coef()%x==0)
            {
                p=(*this).horner(x,bun);
                //if(p.m[p.nr_monoame-1].get_grad()==0 && p.m[p.nr_monoame-1].get_coef()==0)
                if(bun==1)
                {
                    if(x>0)
                        out<<"(x-"<<x<<")*(";
                    else
                        out<<"(x+"<<-x<<")*(";
                    for(i=0;i<p.nr_monoame;i++)
                    {
                        grad=p.m[i].get_grad();
                        coef=p.m[i].get_coef();
                        if(i!=p.nr_monoame-1)
                            out<<coef<<"x^"<<grad<<"+";
                        else
                            out<<coef<<"x^"<<grad;
                    }
                    out<<")"<<endl;
                    return ;
                }
            }
        }
    }

}
void Polinom::initializare()
{
    int i,grad;
    float coef;
    for(i=0;i<nr_monoame;i++)
    {
        cout<<"Grad:";
        cin>>grad;
        cout<<"Coeficient:";
        cin>>coef;
        m[i].set_grad(grad);
        m[i].set_coef(coef);
    }
}
istream& operator>>(istream& in, Polinom &p)
{
    int i,grad;
    float coef;
    for(i=0;i<p.nr_monoame;i++)
    {
        in>>grad>>coef;
        p.m[i].set_grad(grad);
        p.m[i].set_coef(coef);
    }
    return in;
}
ostream& operator<<(ostream& out, Polinom &p)
{
    p.afisare(out);
    return out;
};

int Polinom::crit_eisenstein()
{
    if(m[nr_monoame-1].get_grad()>0)
        return 0;
    int i,x,p;
    for(x=2;x<=m[nr_monoame-1].get_coef();x++)
    {
        p=1;
        for(i=1;i<nr_monoame;i++)
            if(m[i].get_coef()%x!=0)
                p=0;
        if(p==1 && prim(x)==1 && (m[0].get_coef())%x!=0 && (m[nr_monoame-1].get_coef())%(x*x)!=0)
            return 1;
    }
    return 0;
}

int main()
{
    int opt,nr_polinoame,nr_monoame,red,index;
    Polinom* VP[50];
    cout<<"Cate polinoame doriti sa introduceti?"<<endl;
    cin>>nr_polinoame;
    cout<<"Pentru fiecare polinom veti alege daca doriti sa fie reductibil sau ireductibil"<<endl;
    for(int i=0;i<nr_polinoame;i++)
    {
        cout<<"Polinomul cu numarul de ordine "<<i<<" este reductibil(1) sau ireductibil(2)?"<<endl;
        cin>>red;
        if(red==1)
        {
            cout<<"Cate monoame doriti sa contina?"<<endl;
            cin>>nr_monoame;
            VP[i]=new Polinom_reductibil(nr_monoame);
            VP[i]->initializare();
        }
        if(red==2)
        {
            cout<<"Cate monoame doriti sa contina?"<<endl;
            cin>>nr_monoame;
            VP[i]=new Polinom_ireductibil(nr_monoame);
            VP[i]->initializare();
        }
    }
    while(1)
    {
        cout<<"Meniu:"<<endl;
        cout<<"1.Afisare polinom"<<endl;
        cout<<"2.Testare criteriu eisenstein pe polinom"<<endl;
        cin>>opt;
        if(opt==1)
        {
            cout<<"Pe care dintre cele "<<nr_polinoame<<" il alegeti?(primul este 0, ultimul este "<<nr_polinoame-1<<endl;
            cin>>index;
            cout<<(*VP[index])<<endl;
        }
        if(opt==2)
        {
            cout<<"Pe care dintre cele "<<nr_polinoame<<" il alegeti?(primul este 0, ultimul este "<<nr_polinoame-1<<endl;
            cin>>index;
            cout<<VP[index]->crit_eisenstein()<<endl;
        }
    }
    return 0;
}
