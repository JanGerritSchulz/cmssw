#ifndef RecoVertex_Vega_interface_alpaka_PortableOneToManyAssoc_h
#define RecoVertex_Vega_interface_alpaka_PortableOneToManyAssoc_h

// Owning, runtime-sized wrapper around a cms::alpakatools::OneToManyAssoc* type.
//
// The association classes (OneToManyAssocSequential / OneToManyAssocRandomAccess)
// instantiated with kDynamicSize are device-side headers: they contain pointers
// and counters, but own no memory. Someone has to allocate
//   - the header itself   (0-dimensional device buffer, so kernels can mutate it),
//   - the offsets array   (nOnes + 1 counters),
//   - the content array   (capacity of "many" entries),
// and wire the header to the two arrays. This class does that once, in one place.

#include <cstdint>
#include <utility>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/OneToManyAssoc.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  // TAssoc: e.g.
  //   cms::alpakatools::OneToManyAssocSequential<uint32_t, kDynamicSize, kDynamicSize>
  //   cms::alpakatools::OneToManyAssocRandomAccess<uint32_t, kDynamicSize, kDynamicSize>
  template <typename TAssoc>
  class PortableOneToManyAssoc {
  public:
    using Assoc = TAssoc;
    using Storage = typename Assoc::value_type;  // element type of the content array
    using Offsets = typename Assoc::Counter;     // element type of the offsets array
    using View = typename Assoc::View;           // non-owning bundle of pointers and sizes

    // There is no default constructor (device buffers are not default-constructible
    // in general). If you need deferred initialization, hold a
    // std::optional<PortableOneToManyAssoc<...>> and emplace() it once the sizes are known.
    PortableOneToManyAssoc() = delete;

    // nOnes:    number of "ones" (e.g. number of tracks)
    // capacity: maximum total number of "many" entries (e.g. total number of hit-to-track links)
    template <typename TQueue>
    PortableOneToManyAssoc(TQueue& queue, uint32_t nOnes, uint32_t capacity)
        : nOnes_{nOnes},
          capacity_{capacity},
          assoc_{cms::alpakatools::make_device_buffer<Assoc>(queue)},
          offsets_{cms::alpakatools::make_device_buffer<Offsets[]>(queue, nOnes + 1)},
          storage_{cms::alpakatools::make_device_buffer<Storage[]>(queue, capacity)},
          // VERIFY: member order of View. This mirrors the existing CA code:
          //   {header*, offsets*, storage*, nOnes + 1, capacity}
          view_{assoc_.data(), offsets_.data(), storage_.data(), nOnes + 1u, capacity} {
      // Zero the offsets and let the header point at the two arrays.
      Assoc::template launchZero<Acc1D>(view_, queue);
    }

    // Owns device buffers and the view points into them: not copyable.
    PortableOneToManyAssoc(PortableOneToManyAssoc const&) = delete;
    PortableOneToManyAssoc& operator=(PortableOneToManyAssoc const&) = delete;

    // Move: the device buffers are shared_ptr-backed handles, so moving them
    // transfers ownership without moving the underlying device memory. But the
    // view holds raw pointers into those buffers, so it must be REBUILT from
    // the moved buffers - a defaulted move would copy the old view verbatim,
    // leaving it pointing at the moved-from object's buffers (correct only by
    // luck of the buffers being heap handles). Rebuilding restores the
    // invariant "view_ points at my own buffers" by construction.
    PortableOneToManyAssoc(PortableOneToManyAssoc&& other) noexcept
        : nOnes_{other.nOnes_},
          capacity_{other.capacity_},
          assoc_{std::move(other.assoc_)},
          offsets_{std::move(other.offsets_)},
          storage_{std::move(other.storage_)},
          view_{assoc_.data(), offsets_.data(), storage_.data(), nOnes_ + 1u, capacity_} {}

    PortableOneToManyAssoc& operator=(PortableOneToManyAssoc&& other) noexcept {
      if (this != &other) {
        nOnes_ = other.nOnes_;
        capacity_ = other.capacity_;
        assoc_ = std::move(other.assoc_);
        offsets_ = std::move(other.offsets_);
        storage_ = std::move(other.storage_);
        view_ = View{assoc_.data(), offsets_.data(), storage_.data(), nOnes_ + 1u, capacity_};
      }
      return *this;
    }

    // Device pointer to the header: pass this to kernels that fill or read the association.
    Assoc* data() { return assoc_.data(); }
    Assoc const* data() const { return assoc_.data(); }

    // Non-owning view, e.g. for launchFinalize / launchZero helpers taking a View.
    // Valid only as long as this object is alive.
    View const& view() const { return view_; }
    View& view() { return view_; }

    uint32_t nOnes() const { return nOnes_; }
    uint32_t capacity() const { return capacity_; }

  private:
    uint32_t nOnes_;
    uint32_t capacity_;
    // Declaration order matters: the buffers must be constructed before the view uses them.
    cms::alpakatools::device_buffer<Device, Assoc> assoc_;
    cms::alpakatools::device_buffer<Device, Offsets[]> offsets_;
    cms::alpakatools::device_buffer<Device, Storage[]> storage_;
    View view_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif  // RecoVertex_Vega_interface_alpaka_PortableOneToManyAssoc_h
