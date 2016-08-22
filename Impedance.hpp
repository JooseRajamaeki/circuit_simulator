/*
 * Impedance.hpp
 *
 *  Created on: Aug 2, 2016
 *      Author: joose
 */

#ifndef IMPEDANCE_HPP_
#define IMPEDANCE_HPP_

#include <cmath>
#include <complex>

template<typename Scalar>
class Impedance{


	//Frequency dependent
	std::complex<Scalar> impedance_;

	//Non frequency dependent components
	Scalar resistance_;
	Scalar capacitance_;
	Scalar inductance_;

public:


	Impedance(){
		//Open circuit
		impedance_ = (std::complex<Scalar>)0;
		resistance_ = 0;
		capacitance_ = 0;
		inductance_ = 0;
	}

	Impedance(Scalar resistance, Scalar inductance, Scalar capacitance){
		//Open circuit
		impedance_ = (std::complex<Scalar>)0;
		resistance_ = resistance;
		capacitance_ = capacitance;
		inductance_ = inductance;
	}

	Impedance(const Impedance& other){
		impedance_ = other.impedance_;
		resistance_ = other.resistance_;
		capacitance_ = other.capacitance_;
		inductance_ = other.inductance_;
	}

	Impedance operator=(const Impedance& other){
		impedance_ = other.impedance_;
		resistance_ = other.resistance_;
		capacitance_ = other.capacitance_;
		inductance_ = other.inductance_;
		return *this;
	}

	virtual ~Impedance(){

	}

	std::complex<Scalar> get_impedance(void){
		return impedance_;
	}


	void set_resistance(const Scalar& resistance){
		resistance_ = resistance;
	}

	void set_capacitance(const Scalar& capacitance){
		capacitance_ = capacitance;
	}

	void set_inductance(const Scalar& inductance){
		inductance_ = inductance;
	}

	Scalar get_resistance(void){
		return resistance_;
	}

	Scalar get_capacitance(void){
		return capacitance_;
	}

	Scalar get_inductance(void){
		return inductance_;
	}

	virtual void compute_impedance(const Scalar& frequency){
		std::complex<Scalar> imaginary_unit = -1;
		imaginary_unit = std::sqrt(imaginary_unit);
		std::complex<Scalar> omega_i = imaginary_unit * (std::complex<Scalar>)(2 * M_PI * frequency);

		impedance_ = resistance_;
		impedance_ += omega_i * (std::complex<Scalar>)inductance_;
		if (capacitance_ > (Scalar)0){
			impedance_ += (std::complex<Scalar>)1/((std::complex<Scalar>)capacitance_*omega_i);
		}

	}


};





#endif /* IMPEDANCE_HPP_ */
