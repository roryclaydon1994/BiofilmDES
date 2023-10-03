#ifndef INTERFACE_BACTERIA_HPP
#define INTERFACE_BACTERIA_HPP

// Standard modules
#include <iostream>

#include <unordered_map>

class IBacterium
{
public:

  virtual std::string getMyType() const=0;

  virtual Vec3 getPos()       const=0;
  virtual void setPos(double,double,double z=0.0)=0;
  virtual void setAngles(double,double=0.5*constants::pi)=0;

  virtual Vec3 getLoggedPos() const=0;
  virtual void setLoggedPos()=0;
  virtual double getRadius()  const=0;
  virtual double getLength()  const=0;
  virtual double getModE()    const=0;

#ifdef CHAINING
  virtual IBacterium* getUpperLink() const=0;
  virtual IBacterium* getLowerLink() const=0;
  virtual void setUpperLink(IBacterium*)=0;
  virtual void setLowerLink(IBacterium*)=0;
#endif

  virtual void move(double)=0;
  virtual void grow(double)=0;
  virtual void reset()=0;
  virtual void divide(std::vector<IBacterium*>& cell_list)=0;
  virtual bool signalDivide()=0;
  virtual Vec3 getOrientation() const=0;
  virtual void getMyEndVecs(Vec3&,Vec3&) const=0;
  virtual double getEffectiveR() const=0;
  virtual std::vector<IBacterium*>& getNeighbourList()=0;
  virtual IBacterium* getHeadLink()=0;
  virtual void setHeadLink(IBacterium*)=0;
  virtual uint& getGridCell()=0;
  virtual uint getID() const=0;

  virtual Vec3& getForce()=0;
  virtual void addForce(Vec3)=0;
  virtual Vec3& getTorque()=0;
  virtual void addTorque(Vec3)=0;
  virtual Vec3& getVel()=0;
  virtual void setVel()=0;
  virtual Vec3& getAngVel()=0;
  virtual void setAngVel()=0;

  virtual void printToFile(std::ostream &out)=0;

  virtual ~IBacterium () {}

};

inline double getEffectiveQ(const double a, const double b)
{
  return ( a * b ) / ( a + b );
}

#endif //End fileguard
