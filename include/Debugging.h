// C++ dependencies
#include <functional>
#include <iostream>

/// @brief Returns a callback function that prints to std::cout when the value changes.
/// @tparam T The type of the value being watched.
/// @param name The name of the value being watched.
/// @return A callback function that prints to std::cout when the value changes.
template<typename T>
auto coutCallback(const std::string& name)
{
    return [name](const T& oldVal, const T& newVal)
    {
        std::cout << name << " changed from " << oldVal << " to " << newVal << "\n";
    };
}

/// @brief A class that watches a value and calls a callback function when the value changes.
/// @tparam T The type of the value to watch.
template<typename T>
class Watched
{
public:
    using Callback = std::function<void(const T& oldVal, const T& newVal)>;
    
    /// @brief Default constructor.
    Watched() = default;

    /// @brief Default destructor.
    ~Watched()
    {
        if (m_callback) m_callback(m_value, m_value);
        std::cout << "Watched object destroyed\n";
    }

    /// @brief Constructs a Watched object with the given value.
    /// @param value The value to watch.
    Watched(const T& value) : m_value(value) {}

    /// @brief Sets the value being watched.
    /// @param value The value to set.
    /// @return A reference to the watched value.
    Watched& operator=(const T& value) 
    {
        if (m_value != value)
        {
            T old = m_value;
            m_value = value;
            if (m_callback) m_callback(old, m_value);
        }
        return *this;
    }

    /// @brief Returns the value being watched.
    /// @return The value being watched.
    operator const T&() const { return m_value; }

    /// @brief Returns the value being watched.
    /// @return The value being watched.
    const T& get() const { return m_value; }

    /// @brief Sets the callback function to be called when the value changes.
    /// @param cb The callback function.
    void setCallback(Callback cb) { m_callback = std::move(cb); }

private:
    T m_value;
    Callback m_callback;
};