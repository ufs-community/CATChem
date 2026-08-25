#pragma once

namespace catchem {

#ifdef USE_REAL8
    using fp = double;
#else
    using fp = float;
#endif

} // namespace catchem
