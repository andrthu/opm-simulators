/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2019 Dr. Blatt - HPC-Simulation-Software & Services.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
#define OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED


#include <opm/simulators/linalg/twolevelmethodcpr.hh>


namespace Opm
{

template <class FineOperator, class CoarseOperator, class Communication, bool transpose = false>
class PressureTransferPolicy : public Dune::Amg::LevelTransferPolicyCpr<FineOperator, CoarseOperator>
{
public:
    typedef Dune::Amg::LevelTransferPolicyCpr<FineOperator, CoarseOperator> ParentType;
    typedef Communication ParallelInformation;
    typedef typename FineOperator::domain_type FineVectorType;

public:
    PressureTransferPolicy(const Communication& comm, const FineVectorType& weights, int pressure_var_index)
        : communication_(&const_cast<Communication&>(comm))
        , weights_(weights)
        , pressure_var_index_(pressure_var_index)
    {
    }

    virtual void createCoarseLevelSystem(const FineOperator& fineOperator) override
    {
        using CoarseMatrix = typename CoarseOperator::matrix_type;
        const auto& fineLevelMatrix = fineOperator.getmat();
	setIS(fineOperator);
	interiorSize_ = fineLevelMatrix.N();
	if (communication_->communicator().rank() == 0) {std::cout<< interiorSize_<<" "<< fineLevelMatrix.N() << std::endl;}
        coarseLevelMatrix_.reset(new CoarseMatrix(fineLevelMatrix.N(), fineLevelMatrix.M(), CoarseMatrix::row_wise));
        auto createIter = coarseLevelMatrix_->createbegin();

        for (const auto& row : fineLevelMatrix) {
            for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
                createIter.insert(col.index());
            }
            ++createIter;
        }

	//calculateCoarseGhostEntries(fineLevelMatrix.N());
        calculateCoarseEntries(fineOperator);
        coarseLevelCommunication_.reset(communication_, [](Communication*) {});

        this->lhs_.resize(this->coarseLevelMatrix_->M());
        this->rhs_.resize(this->coarseLevelMatrix_->N());
        using OperatorArgs = typename Dune::Amg::ConstructionTraits<CoarseOperator>::Arguments;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
        OperatorArgs oargs(coarseLevelMatrix_, *coarseLevelCommunication_);
        this->operator_ = Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs);
#else
        OperatorArgs oargs(*coarseLevelMatrix_, *coarseLevelCommunication_);
        this->operator_.reset(Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs));
#endif
    }

    virtual void calculateCoarseEntries(const FineOperator& fineOperator) override
    {
        const auto& fineMatrix = fineOperator.getmat();
        *coarseLevelMatrix_ = 0;
        auto rowCoarse = coarseLevelMatrix_->begin();
        //for (auto row = fineMatrix.begin(), rowEnd = fineMatrix.end(); row != rowEnd; ++row, ++rowCoarse) {
	for (auto row = fineMatrix.begin(); row.index() < interiorSize_; ++row, ++rowCoarse) {
            assert(row.index() == rowCoarse.index());
            auto entryCoarse = rowCoarse->begin();
            for (auto entry = row->begin(), entryEnd = row->end(); entry != entryEnd; ++entry, ++entryCoarse) {
                assert(entry.index() == entryCoarse.index());
                double matrix_el = 0;
                if (transpose) {
                    const auto& bw = weights_[entry.index()];
                    for (size_t i = 0; i < bw.size(); ++i) {
                        matrix_el += (*entry)[pressure_var_index_][i] * bw[i];
                    }
                } else {
                    const auto& bw = weights_[row.index()];
                    for (size_t i = 0; i < bw.size(); ++i) {
                        matrix_el += (*entry)[i][pressure_var_index_] * bw[i];
                    }
                }
                (*entryCoarse) = matrix_el;
            }
        }
        //assert(rowCoarse == coarseLevelMatrix_->end());
	calculateCoarseGhostEntries(fineMatrix.N());
    }

    void calculateCoarseGhostEntries(size_t numRows)
    {
	for (size_t row = interiorSize_; row < numRows; ++row) {
	    (*coarseLevelMatrix_)[row][row] = 1.0;
	}
    }

    virtual void moveToCoarseLevel(const typename ParentType::FineRangeType& fine) override
    {
        // Set coarse vector to zero
        this->rhs_ = 0;

        auto begin = fine.begin();

        for (auto block = begin; block.index() < interiorSize_; ++block) {
            const auto& bw = weights_[block.index()];
            double rhs_el = 0.0;
            if (transpose) {
                rhs_el = (*block)[pressure_var_index_];
            } else {
                for (size_t i = 0; i < block->size(); ++i) {
                    rhs_el += (*block)[i] * bw[i];
                }
            }
            this->rhs_[block - begin] = rhs_el;
        }

        this->lhs_ = 0;
    }

    virtual void moveToFineLevel(typename ParentType::FineDomainType& fine) override
    {
        //auto end = fine.end();
	auto begin = fine.begin();

        for (auto block = begin; block.index() < interiorSize_; ++block) {
            if (transpose) {
                const auto& bw = weights_[block.index()];
                for (size_t i = 0; i < block->size(); ++i) {
                    (*block)[i] = this->lhs_[block - begin] * bw[i];
                }
            } else {
                (*block)[pressure_var_index_] = this->lhs_[block - begin];
            }
        }
    }

    virtual PressureTransferPolicy* clone() const override
    {
        return new PressureTransferPolicy(*this);
    }

    const Communication& getCoarseLevelCommunication() const
    {
        return *coarseLevelCommunication_;
    }

    std::size_t getPressureIndex() const
    {
        return pressure_var_index_;
    }
private:

    template<class O>
    void setIS(const O& ad) {interiorSize_ = ad.getmat().N();}

    void setIS(const GhostLastMatrixAdapter<typename FineOperator::matrix_type,FineVectorType,FineVectorType,Communication>& ad){interiorSize_ = ad.getInteriorSize();}
    Communication* communication_;
    const FineVectorType& weights_;
    const std::size_t pressure_var_index_;
    std::shared_ptr<Communication> coarseLevelCommunication_;
    std::shared_ptr<typename CoarseOperator::matrix_type> coarseLevelMatrix_;
    std::size_t interiorSize_;
};

} // namespace Opm

#endif // OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
