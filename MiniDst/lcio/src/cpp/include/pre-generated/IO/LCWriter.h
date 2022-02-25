// -*- C++ -*-
// AID-GENERATED
// =========================================================================
// This class was generated by AID - Abstract Interface Definition          
// DO NOT MODIFY, but use the org.freehep.aid.Aid utility to regenerate it. 
// =========================================================================
#ifndef IO_LCWRITER_H
#define IO_LCWRITER_H 1

#include <string>

#include "EVENT/LCEvent.h"
#include "EVENT/LCRunHeader.h"
#include "Exceptions.h"

namespace IO {

/**Interface for writing data with LCIO. Uses interfaces 
 * from EVENT/hep.lcio.event.
 * Use LCFactory to instantiate a corresponding LCWriter object
 * for the output format at hand (SIO only, so far).
 *
 * @see LCFactory
 * @author gaede
 * @version Mar 4, 2003
 */
class LCWriter {

public: 
    /// Destructor.
    virtual ~LCWriter() { /* nop */; }

    /** Opens a file for writing. If file with given name exists, 
     * an exception is thrown. Use append or new mode instead.
     *
     *@throws IOException
     */
    virtual void open(const std::string & filename) = 0;

    /** Opens a file for writing.
     * Possible write modes are: LCIO::WRITE_NEW
     * (existing files are replaced) and LCIO::WRITE_APPEND. 
     *
     *@throws IOException
     */
    virtual void open(const std::string & filename, int writeMode) = 0;

    /** Set the compression level - needs to be called before open() otherwise
     *  call will have no effect. If not called the Writer will use default compression.<br>
     *  Valid compression levels are:
     *  <ul>
     *    <li> level <  0 : default compression </li>
     *    <li> level == 0 : no compression</li>
     *    <li> level >  0 : compression level (typically 1 (fastest) - 9 (best compression))
     *    </li>
     *  </ul>
     *  
     *  Status:  (v01-09)<br>
     *  C++: experimental code - don't use for production<br>
     *  Java: not implemented
     * 
     *@param level compression level
     */
    virtual void setCompressionLevel(int level) = 0;

    /** Writes the given run header to file.
     *
     *@throws IOException
     */
    virtual void writeRunHeader(const EVENT::LCRunHeader * hdr)  = 0;

    /** Writes the given event to file.
     *
     *@throws IOException
     */
    virtual void writeEvent(const EVENT::LCEvent * evt) = 0;

    /** Closes the output file/stream.
     *
     *@throws IOException
     */
    virtual void close() = 0;

    /** Flushes the output file/stream.
     *
     *@throws IOException
     */
    virtual void flush()  = 0;
}; // class
} // namespace IO
#endif /* ifndef IO_LCWRITER_H */
