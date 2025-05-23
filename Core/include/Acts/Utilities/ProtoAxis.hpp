// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <memory>
#include <vector>

namespace Acts {
/// This class is a pure run-time placeholder for the axis definition
///
/// The IAxis allows via the visitor pattern to access the actual axis type
/// which helps to create grid creation code by the compiler as done
/// in the makeGrid helper functions.
///
/// In addition to a simple axis definitions, it holds also a description
/// of the axis direction.
class ProtoAxis {
 public:
  /// Convenience constructors - for variable binning
  ///
  /// @param abType the axis boundary type
  /// @param edges the bin edges (variable binning)
  ProtoAxis(Acts::AxisBoundaryType abType, const std::vector<double>& edges);

  /// Convenience constructors - for equidistant binning
  ///
  /// @param abType the axis boundary type
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  /// @param nbins the number of bins
  ProtoAxis(AxisBoundaryType abType, double minE, double maxE,
            std::size_t nbins);

  /// Placeholder constructor for auto-range binning
  ///
  /// @param abType the axis boundary type
  /// @param nbins the number of bins
  ///
  /// @note that auto-range is only supported for equidistant binning
  ProtoAxis(AxisBoundaryType abType, std::size_t nbins);

  /// Custom copy constructor
  /// @param other is the right hand side ProtoAxis
  ProtoAxis(const ProtoAxis& other);

  /// Custom assignment operator
  /// @param other is the right hand side ProtoAxis
  ProtoAxis& operator=(const ProtoAxis& other);

  ProtoAxis(ProtoAxis&&) = default;
  ProtoAxis& operator=(ProtoAxis&&) = default;

  ~ProtoAxis() = default;

  /// @brief Return the IAxis representation
  ///
  /// @return @c AxisType of this axis
  const IAxis& getAxis() const;

  /// Set the range, in case of autorange, this
  /// will toggle the autorange flag to false
  ///
  /// @throws an exception if minE > maxE
  ///
  /// @note In case of variable binning, it will clip the bins outside
  /// the new range off, i.e. it will potentially change the number
  /// of bins.
  ///
  /// @note In case of equidistant binning, it will adapt the bin size,
  /// and NOT change the number of bins.
  ///
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  void setRange(double minE, double maxE);

  /// @brief check if this is an auto-range binning
  bool isAutorange() const;

  /// Dump into a string
  /// @return the string representation
  std::string toString() const;

 private:
  /// @brief Hidden friend ostream operator
  /// @param os the output stream
  /// @param a the proto axis
  /// @return the streamed into operator
  friend std::ostream& operator<<(std::ostream& os, const ProtoAxis& a) {
    os << a.toString();
    return os;
  }

  /// Dispatch to the correct stream operator
  /// @param os output stream
  void toStream(std::ostream& os) const;

 protected:
  /// The axis representation
  std::unique_ptr<IAxis> m_axis = nullptr;

  /// Indicate if this is a place holder auto-range binning
  bool m_autorange = false;
};

/// @brief Helper method to create a 1D grid from a single proto axis
///
/// @tparam payload_t the grid payloat type
///
/// @param a the proto axis
///
/// @return an IGrid unique ptr and hence transfers ownership
template <typename payload_t>
std::unique_ptr<IGrid> makeGrid(const ProtoAxis& a) {
  if (a.isAutorange()) {
    throw std::invalid_argument(
        "ProtoAxis::makeGrid: Auto-range of the proto axis is not (yet) "
        "resolved, call setRange() first.");
  }

  return a.getAxis().visit(
      [&]<typename AxisTypeA>(const AxisTypeA& axis) -> std::unique_ptr<IGrid> {
        using GridType = Grid<payload_t, AxisTypeA>;
        return std::make_unique<GridType>(axis);
      });
}

/// @brief Helper method to create a 2D grid from a two proto axes
///
/// @tparam payload_t the grid payloat type
///
/// @param a the first proto axis
/// @param b the second proto axis
///
/// @return an IGrid unique ptr and hence transfers ownership
template <typename payload_t>
std::unique_ptr<IGrid> makeGrid(const ProtoAxis& a, const ProtoAxis& b) {
  if (a.isAutorange() || b.isAutorange()) {
    throw std::invalid_argument(
        "ProtoAxis::makeGrid: Auto-range of the proto axis is not (yet) "
        "resolved, call setRange() first.");
  }

  return a.getAxis().visit([&]<typename AxisTypeA>(const AxisTypeA& axisA)
                               -> std::unique_ptr<IGrid> {
    return b.getAxis().visit([&]<typename AxisTypeB>(const AxisTypeB& axisB)
                                 -> std::unique_ptr<IGrid> {
      using GridType = Grid<payload_t, AxisTypeA, AxisTypeB>;
      return std::make_unique<GridType>(axisA, axisB);
    });
  });
}

/// A Directed proto axis
struct DirectedProtoAxis : public ProtoAxis {
 public:
  /// Convenience constructors - for variable binning
  ///
  /// @param axisDir the axis direction
  /// @param abType the axis boundary type
  /// @param edges the bin edges (variable binning)
  DirectedProtoAxis(AxisDirection axisDir, AxisBoundaryType abType,
                    const std::vector<double>& edges);

  /// Convenience constructors - for equidistant binning
  ///
  /// @param axisDir the axis direction
  /// @param abType the axis boundary type
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  /// @param nbins the number of bins
  DirectedProtoAxis(AxisDirection axisDir, AxisBoundaryType abType, double minE,
                    double maxE, std::size_t nbins);

  /// Placeholder constructor for auto-range binning
  ///
  /// @param axisDir the axis direction
  /// @param abType the axis boundary type
  /// @param nbins the number of bins
  ///
  /// @note that auto-range is only supported for equidistant binning
  DirectedProtoAxis(AxisDirection axisDir, AxisBoundaryType abType,
                    std::size_t nbins);

  virtual ~DirectedProtoAxis() = default;

  /// Access to the AxisDirection
  /// @return the axis direction
  AxisDirection getAxisDirection() const;

  /// Dump into a string
  /// @return the string representation
  std::string toString() const;

 private:
  /// @brief Hidden friend ostream operator
  /// @param os the output stream
  /// @param a the proto axis
  /// @return the streamed into operator
  friend std::ostream& operator<<(std::ostream& os,
                                  const DirectedProtoAxis& a) {
    os << a.toString();
    return os;
  }

  /// Dispatch to the correct stream operator
  /// @param os output stream
  void toStream(std::ostream& os) const;

  AxisDirection m_direction;
};

std::ostream& operator<<(std::ostream& os, const std::vector<ProtoAxis>& a);

std::ostream& operator<<(std::ostream& os,
                         const std::vector<DirectedProtoAxis>& a);

}  // namespace Acts
