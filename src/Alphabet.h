#ifndef ALPHABET_H
#define ALPHABET_H

#include "main.h"
#include "debug.h"
#include <string>
#include <assert.h>

class Alphabet {
protected:
	Alphabet() {
		data = -1;
	};
	char data;
public:
	enum {
		DIM = 0
	};

	static const Alphabet GAP;
	static const Alphabet X;
	static const Alphabet stripStart;
	static const Alphabet stripEnd;
	int value() const;
	std::string asString();
	char asChar() const;
	bool isGap() const;
	bool isValid() const;
	bool isUnknown() const;
//	operator char() const;
//	operator int() const;
	bool operator==(const Alphabet &other) const;
	bool operator!=(const Alphabet &other) const;
	bool operator<(const Alphabet &other) const;
};

class AA : Alphabet {
public:
	explicit AA(char c);
	explicit AA(int i);
	AA() {};

	enum {
		DIM = 20
	};

	static const AA GAP;
	static const AA X;
	static const AA stripStart;
	static const AA stripEnd;
	int value() const;
	std::string asString() const;
	char asChar() const;
	bool isGap() const { return this->data == this->GAP.data; }
	bool isValid() const { return this->value() >= 0 && this->value() < this->DIM; }
	bool isUnknown() const { return this->data == this->X.data; }
//	operator char() const { return this->asChar(); }
//	operator int() const { return this->value(); }
	bool operator==(const AA &other) const { return this->data == other.data; }
	bool operator!=(const AA &other) const { return this->data != other.data; }
	bool operator<(const AA &other) const { return this->data < other.data; }
};

class Codon : Alphabet {
public:
	explicit Codon(char c1, char c2, char c3);
	explicit Codon(int i);
	Codon() {};

	enum {
		DIM = 61
	};

	static const Codon GAP;
	static const Codon X;
	static const Codon stripStart;
	static const Codon stripEnd;
	int value() const;
	std::string asString() const;
	char asChar() const;
	operator AA() const { return AA(this->asChar()); }
	bool isGap() const { return this->data == this->GAP.data; }
	bool isValid() const { return this->value() >= 0 && this->value() < this->DIM; }
	bool isUnknown() const { return this->data == this->X.data; }
//	operator char() const { return this->asChar(); }
//	operator int() const { return this->value(); }
	bool operator==(const Codon &other) const { return this->data == other.data; }
	bool operator!=(const Codon &other) const { return this->data != other.data; }
	bool operator<(const Codon &other) const { return this->data < other.data; }
};

class DNA : Alphabet {
public:
	explicit DNA(char c);
	explicit DNA(int i);
	DNA() {};

	enum {
		DIM = 4
	};

	static const DNA GAP;
	static const DNA X;
	static const DNA stripStart;
	static const DNA stripEnd;
	int value() const;
	std::string asString() const;
	char asChar() const;
	bool isGap() const { return this->data == this->GAP.data; }
	bool isValid() const { return this->value() >= 0 && this->value() < this->DIM; }
	bool isUnknown() const { return this->data == this->X.data; }
//	operator char() const { return this->asChar(); }
//	operator int() const { return this->value(); }
	bool operator==(const DNA &other) const { return this->data == other.data; }
	bool operator!=(const DNA &other) const { return this->data != other.data; }
	bool operator<(const DNA &other) const { return this->data < other.data; }
};


#define sequence_t std::basic_string

template<class T>
sequence_t<T> sequenceFromString (const std::string &str) {
	sequence_t<T> seq;
	seq.reserve(str.length());

	for(index_t i=0; i<str.length(); ++i) {
		T c = T(str[i]);
		if(c == T::GAP) {
			error("No support for gapped sequences (yet)");
		}
		seq += c;
	}

	return seq;
}

template<> sequence_t<Codon> sequenceFromString<Codon>(const std::string &str);

template<class T>
std::string stringFromSequence(const sequence_t<T> &seq) {
	std::string str;
	str.reserve(seq.length());

	for(index_t i=0; i<seq.length(); ++i) {
		str += seq[i].asChar();
	}

	return str;
}

template<class T>
std::string stringFromSequence(const sequence_t<T> &seq, const std::string &orig) {
	std::string str;
	str.reserve(orig.length());

	index_t j=0;
	for(index_t i=0; i<seq.length(); ++i) {
		if(seq[i].isGap()) {
			str += seq[i].asChar();
		} else {
			str += orig[j++];
		}
	}

	assert(j == orig.length());

	return str;
}

template<> std::string stringFromSequence<Codon>(const sequence_t<Codon> &seq);
template<> std::string stringFromSequence<Codon>(const sequence_t<Codon> &seq, const std::string &orig);

sequence_t<AA> translateCodons(const sequence_t<Codon> &seq);

#endif /* ALPHABET_H */
